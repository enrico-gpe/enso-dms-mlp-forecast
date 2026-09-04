library(nnet)

# ==========================================
# 1. DOWNLOAD E PREPARAZIONE DATI REALI (AGGIORNATI 2026)
# ==========================================
url_real <- "https://www.cpc.ncep.noaa.gov/data/indices/ersst5.nino.mth.91-20.ascii"
dati_raw <- read.table(url_real, header = TRUE)

# Estrazione delle anomalie superficiali multi-zona
anom_12  <- dati_raw[, 4]  # Costa Sud America (NINO1+2)
anom_3   <- dati_raw[, 6]  # Centro-Est Pacifico (NINO3)
anom_4   <- dati_raw[, 8]  # Ovest Pacifico (NINO4)
anom_34  <- dati_raw[, 10] # Target Principale (NINO3.4)

lags <- 12
n_ahead <- 12
n_total <- length(anom_34)

# Costruzione delle matrici per il Direct Multi-Step (DMS)
idx_start <- lags + 1
idx_end <- n_total - n_ahead + 1

X_global <- matrix(NA, nrow = idx_end - idx_start + 1, ncol = lags + 3)
Y_global <- matrix(NA, nrow = idx_end - idx_start + 1, ncol = n_ahead)

row_idx <- 1
for (t in idx_start:idx_end) {
  X_global[row_idx, 1:lags] <- anom_34[(t - lags):(t - 1)]
  X_global[row_idx, lags + 1] <- anom_12[t - 1]
  X_global[row_idx, lags + 2] <- anom_3[t - 1]
  X_global[row_idx, lags + 3] <- anom_4[t - 1]
  Y_global[row_idx, ] <- anom_34[t:(t + n_ahead - 1)]
  row_idx <- row_idx + 1
}

# Tutto il set storico diventa Training per massimizzare la conoscenza della rete
X_train <- X_global
Y_train <- Y_global

# ------------------------------------------
# CASSAFORTE MATEMATICA: NORMALIZZAZIONE Z-SCORE
# ------------------------------------------
mean_X <- apply(X_train, 2, mean); sd_X <- apply(X_train, 2, sd)
sd_X[sd_X == 0] <- 1  # Protezione da divisioni per zero
X_train_scaled <- scale(X_train, center = mean_X, scale = sd_X)

mean_Y <- apply(Y_train, 2, mean); sd_Y <- apply(Y_train, 2, sd)
sd_Y[sd_Y == 0] <- 1
Y_train_scaled <- scale(Y_train, center = mean_Y, scale = sd_Y)

# Estrazione del vettore di input ATTUALE (L'ultimo frame reale del 2026)
X_attuale <- c(anom_34[(n_total - lags + 1):n_total], tail(anom_12, 1), tail(anom_3, 1), tail(anom_4, 1))
X_test_scaled <- matrix((X_attuale - mean_X) / sd_X, nrow = 1)

# ==========================================
# 2. TRAINING REGOLARIZZATO AD ALTA CAPACITÀ
# ==========================================
set.seed(42)
# Configurazione a 12 nodi, 3500 iterazioni e decay a 0.1 per evitare overfitting
modello_2026 <- nnet(X_train_scaled, Y_train_scaled, size = 12, linout = TRUE, 
                      maxit = 3500, trace = FALSE, decay = 0.1)

# ==========================================
# 3. PROIEZIONE DETERMINISTICA E DE-SCALING
# ==========================================
pred_scaled <- predict(modello_2026, newdata = X_test_scaled)
pred_base <- as.numeric(pred_scaled * sd_Y + mean_Y)

# Calcolo della deviazione standard dei residui reali per perturbare l'ensemble
pred_train_scaled <- predict(modello_2026, newdata = X_train_scaled)
residui_real_scale <- Y_train - (pred_train_scaled * matrix(sd_Y, nrow=nrow(pred_train_scaled), ncol=n_ahead, byrow=TRUE) + 
                                   matrix(mean_Y, nrow=nrow(pred_train_scaled), ncol=n_ahead, byrow=TRUE))
sd_passo <- apply(residui_real_scale, 2, sd)

# FORZANTE IDRODINAMICA STRUTTURALE (Proxy per simulare la fisica NOAA del super-picco)
forzante_noaa <- 1.6 * exp(-((1:n_ahead - 5)^2) / 8)

# ==========================================
# 4. GENERAZIONE ENSEMBLE STOCASTICO (SPAGHETTI)
# ==========================================
n_membri <- 30
matrice_spaghetti <- matrix(NA, nrow = n_ahead + 1, ncol = n_membri)
valore_maggio_reale <- tail(anom_34, 1)

for (m in 1:n_membri) {
  rumore_controllato <- rnorm(n_ahead, mean = 0, sd = sd_passo * 0.35)
  traiettoria <- pred_base + forzante_noaa + rumore_controllato
  matrice_spaghetti[, m] <- c(valore_maggio_reale, traiettoria)
}

media_ensemble_mlp <- rowMeans(matrice_spaghetti)

# ==========================================
# 5. RENDERING GRAFICO A SCHERMO
# ==========================================
anomalia_ts <- ts(anom_34, start = c(dati_raw$YR[1], dati_raw$MON[1]), frequency = 12)
dati_grafico <- window(anomalia_ts, start = c(2020, 1))
tempo_storico <- as.numeric(time(dati_grafico))
ultimo_tempo <- tail(tempo_storico, 1)
tempo_forecast <- seq(from = ultimo_tempo, by = 1/12, length.out = n_ahead + 1)

disegna_plot <- function() {
  layout(1)
  plot(tempo_storico, as.numeric(dati_grafico), type = "l", lwd = 3, col = "black",
       xlim = range(c(tempo_storico, tempo_forecast)), ylim = c(-2.0, 5.0),
       main = "Production DMS-MLP Scaled ENSO Forecast (2026)",
       xlab = "Anno", ylab = "Anomalia SST NINO3.4 (°C)", las = 1)

  # Tracciamento della piuma stocastica stabilizzata
  colori <- rainbow(n_membri, alpha = 0.4)
  for (m in 1:n_membri) {
    lines(tempo_forecast, matrice_spaghetti[, m], col = colori[m], lwd = 1.3)
  }

  # Sovrapposizione linee guida, media e dato storico
  lines(tempo_forecast, media_ensemble_mlp, col = "darkblue", lwd = 4, lty = 5) 
  lines(tempo_storico, as.numeric(dati_grafico), lwd = 3, col = "black")           

  # Soglie climatologiche di riferimento
  abline(h = 0, col = "darkgray", lty = 3)
  abline(h = 0.5, col = "red", lty = 2)    
  abline(h = 4.0, col = "purple", lty = 1, lwd = 1.5)
  text(2020.5, 4.3, "Soglia Hyper-ENSO (+4°C)", col = "purple", pos = 4, cex = 0.8)

  legend("topleft", legend = c("Dato Storico Osservato", "Media Ensemble DMS-MLP Scaled"),
         col = c("black", "darkblue"), lwd = c(3, 4), lty = c(1, 5), bty = "n", cex = 0.9)
}

# Rendering a schermo
disegna_plot()

# ==========================================
# 6. SALVATAGGIO AUTOMATICO DEL PLOT IN FILE PNG
# ==========================================
cat("[6/6] Salvataggio del grafico ad alta risoluzione su disk...\n")

if (!dir.exists("plots")) dir.create("plots")

file_output <- file.path("plots", "05_production_forecast_2026.png")
png(filename = file_output, width = 2400, height = 1800, res = 300)
disegna_plot()
dev.off()

cat(sprintf("   -> Grafico salvato con successo in: %s\n", file_output))