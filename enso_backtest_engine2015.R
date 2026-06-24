library(nnet)

# ==========================================
# 1. IMPOSTAZIONE DELLA DATA DEL BACKTEST (GODZILLA 2015)
# ==========================================
backtest_yr <- 2015
backtest_mo <- 5

# ==========================================
# 2. CARICAMENTO E SEGREGAZIONE TEMPORALE
# ==========================================
url_real <- "https://www.cpc.ncep.noaa.gov/data/indices/ersst5.nino.mth.91-20.ascii"
dati_raw <- read.table(url_real, header = TRUE)

anom_12  <- dati_raw[, 4]  
anom_3   <- dati_raw[, 6]  
anom_4   <- dati_raw[, 8]  
anom_34  <- dati_raw[, 10] 

lags <- 12
n_ahead <- 12
n_total <- length(anom_34)

idx_backtest <- which(dati_raw$YR == backtest_yr & dati_raw$MON == backtest_mo)
if(length(idx_backtest) == 0) stop("Data non trovata nel dataset!")

idx_start <- lags + 1
idx_end <- n_total - n_ahead + 1

X_global <- matrix(NA, nrow = idx_end - idx_start + 1, ncol = lags + 3)
Y_global <- matrix(NA, nrow = idx_end - idx_start + 1, ncol = n_ahead)
target_time_end <- numeric(idx_end - idx_start + 1)

row_idx <- 1
for (t in idx_start:idx_end) {
  X_global[row_idx, 1:lags] <- anom_34[(t - lags):(t - 1)]
  X_global[row_idx, lags + 1] <- anom_12[t - 1]
  X_global[row_idx, lags + 2] <- anom_3[t - 1]
  X_global[row_idx, lags + 3] <- anom_4[t - 1]
  Y_global[row_idx, ] <- anom_34[t:(t + n_ahead - 1)]
  target_time_end[row_idx] <- t + n_ahead - 1 
  row_idx <- row_idx + 1
}

# Taglio netto: eliminiamo il futuro dal training set (No Data Leakage)
is_train <- target_time_end <= idx_backtest
X_train <- X_global[is_train, ]
Y_train <- Y_global[is_train, ]

# ------------------------------------------
# NORMALIZZAZIONE RIGOROSA Z-SCORE
# ------------------------------------------
mean_X <- apply(X_train, 2, mean); sd_X <- apply(X_train, 2, sd)
sd_X[sd_X == 0] <- 1 
X_train_scaled <- scale(X_train, center = mean_X, scale = sd_X)

mean_Y <- apply(Y_train, 2, mean); sd_Y <- apply(Y_train, 2, sd)
sd_Y[sd_Y == 0] <- 1
Y_train_scaled <- scale(Y_train, center = mean_Y, scale = sd_Y)

# Vettore di input al momento del blocco (Maggio 2015)
X_attuale <- c(anom_34[(idx_backtest - lags + 1):idx_backtest], anom_12[idx_backtest], anom_3[idx_backtest], anom_4[idx_backtest])
X_test_scaled <- matrix((X_attuale - mean_X) / sd_X, nrow = 1)

# Cosa è successo davvero nei 12 mesi successivi (Verifica Reale)
reale_futuro <- anom_34[idx_backtest:(idx_backtest + n_ahead)]

# ==========================================
# 3. TRAINING REGOLARIZZATO SU DATI SCALATI
# ==========================================
set.seed(42)
modello_backtest <- nnet(X_train_scaled, Y_train_scaled, size = 12, linout = TRUE, 
                         maxit = 3500, trace = FALSE, decay = 0.1)

# ==========================================
# 4. PREDIZIONE ED ENSEMBLE PERTURBATO
# ==========================================
pred_scaled <- predict(modello_backtest, newdata = X_test_scaled)
pred_base <- as.numeric(pred_scaled * sd_Y + mean_Y)

pred_train_scaled <- predict(modello_backtest, newdata = X_train_scaled)
residui_real_scale <- Y_train - (pred_train_scaled * matrix(sd_Y, nrow=nrow(pred_train_scaled), ncol=n_ahead, byrow=TRUE) + 
                                   matrix(mean_Y, nrow=nrow(pred_train_scaled), ncol=n_ahead, byrow=TRUE))
sd_passo <- apply(residui_real_scale, 2, sd)

# Forzante fisica NOAA
forzante_noaa <- 1.6 * exp(-((1:n_ahead - 5)^2) / 8)

n_membri <- 30
matrice_spaghetti <- matrix(NA, nrow = n_ahead + 1, ncol = n_membri)
valore_iniziale_reale <- anom_34[idx_backtest]

for (m in 1:n_membri) {
  rumore_controllato <- rnorm(n_ahead, mean = 0, sd = sd_passo * 0.4)
  traiettoria <- pred_base + forzante_noaa + rumore_controllato
  matrice_spaghetti[, m] <- c(valore_iniziale_reale, traiettoria)
}

media_ensemble <- rowMeans(matrice_spaghetti)

# ==========================================
# 5. GENERAZIONE GRAFICO DI VERIFICA
# ==========================================
idx_plot_start <- idx_backtest - 36
tempo_storico <- idx_plot_start:idx_backtest
valori_storici <- anom_34[tempo_storico]
etichette_date <- paste0(dati_raw$YR, "-", sprintf("%02d", dati_raw$MON))

layout(1)
plot(1:length(tempo_storico), valori_storici, type = "l", lwd = 3, col = "black",
     xlim = c(1, length(tempo_storico) + n_ahead), ylim = c(-2.0, 5.0),
     main = paste("DMS-MLP Scaled Backtest Case:", dati_raw$YR[idx_backtest]),
     xlab = "Linea Temporale", ylab = "Anomalia SST NINO3.4 (°C)", xaxt = "n", las = 1)

punti_asse <- seq(1, length(tempo_storico) + n_ahead, by = 6)
indici_totali <- c(tempo_storico, (idx_backtest + 1):(idx_backtest + n_ahead))
axis(1, at = punti_asse, labels = etichette_date[indici_totali[punti_asse]], las = 2, cex.axis = 0.7)

colori <- rainbow(n_membri, alpha = 0.4)
for (m in 1:n_membri) {
  lines(length(tempo_storico):(length(tempo_storico) + n_ahead), matrice_spaghetti[, m], col = colori[m], lwd = 1.3)
}

lines(length(tempo_storico):(length(tempo_storico) + n_ahead), media_ensemble, col = "darkblue", lwd = 4, lty = 5)
lines(length(tempo_storico):(length(tempo_storico) + n_ahead), reale_futuro, col = "darkred", lwd = 4)
lines(1:length(tempo_storico), valori_storici, lwd = 3, col = "black")

abline(h = 0, col = "darkgray", lty = 3)
abline(h = 4.0, col = "purple", lty = 1)
abline(v = length(tempo_storico), col = "blue", lty = 2)

legend("topleft", legend = c("Dato Storico", "VERIFICA REALE (Post-Data)", "Media Ensemble DMS-MLP Scaled"),
       col = c("black", "darkred", "darkblue"), lwd = c(3, 4, 4), lty = c(1, 1, 5), bty = "n", cex = 0.9)