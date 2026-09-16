# Copyright (C) 2026 Enrico Pozzi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

library(nnet)
library(tidyverse)
library(lubridate)

# ==========================================
# 1. PREPARAZIONE DATI
# ==========================================
cat("[1/4] Download e pulizia dati SST in corso da NOAA CPC...\n")

url_real <- "https://www.cpc.ncep.noaa.gov/data/indices/rel_wksst9120.txt"
righe_grezze <- readLines(url_real)
righe_pulite <- trimws(righe_grezze)
idx_dati <- grep("^[0-9]{1,2}[a-zA-Z]{3}[0-9]{4}", righe_pulite, ignore.case = TRUE)

dati_raw <- read.table(text = righe_pulite[idx_dati], header = FALSE, sep = "", fill = TRUE)[, 1:5]
colnames(dati_raw) <- c("Week", "ANOM12", "ANOM3", "ANOM34", "ANOM4")

dati_raw <- dati_raw %>%
  mutate(Week = dmy(Week), across(ANOM12:ANOM4, as.numeric)) %>%
  filter(!is.na(Week)) %>%
  drop_na() %>%
  arrange(Week)

data_target <- ymd("2026-09-02")
ultima_obs <- max(dati_raw$Week)

if (ultima_obs < data_target) {
  wks_missing <- seq(ultima_obs + weeks(1), data_target, by = "1 week")
  n_miss <- length(wks_missing)
  df_miss <- data.frame(Week = wks_missing)
  
  for (col in c("ANOM12", "ANOM3", "ANOM34", "ANOM4")) {
    fit <- lm(y ~ x, data = data.frame(x = 1:4, y = tail(dati_raw[[col]], 4)))
    df_miss[[col]] <- round(as.numeric(predict(fit, newdata = data.frame(x = 4 + (1:n_miss)))), 2)
  }
  dati_raw <- rbind(dati_raw, df_miss)
}

anom_12 <- dati_raw$ANOM12
anom_3  <- dati_raw$ANOM3
anom_4  <- dati_raw$ANOM4
anom_34 <- dati_raw$ANOM34

# ==========================================
# 2. CONFIGURAZIONE MATRICI E SCALING
# ==========================================
lags <- 12
n_ahead <- 18
n_total <- length(anom_34)

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

X_train <- X_global; Y_train <- Y_global

mean_X <- apply(X_train, 2, mean); sd_X <- apply(X_train, 2, sd); sd_X[sd_X == 0] <- 1
X_train_scaled <- scale(X_train, center = mean_X, scale = sd_X)

mean_Y <- apply(Y_train, 2, mean); sd_Y <- apply(Y_train, 2, sd); sd_Y[sd_Y == 0] <- 1
Y_train_scaled <- scale(Y_train, center = mean_Y, scale = sd_Y)

X_attuale <- c(anom_34[(n_total - lags + 1):n_total], tail(anom_12, 1), tail(anom_3, 1), tail(anom_4, 1))

# ==========================================
# 3. GENERAZIONE ENSEMBLE DINAMICO E FRASTAGLIATO
# ==========================================
n_members <- 50
ensemble_matrix_pure <- matrix(NA, nrow = n_members, ncol = n_ahead)

cat(sprintf("[3/4] Addestramento Ensemble Dinamico (%d membri con perturbazione strutturale e stocastica)...\n", n_members))

set.seed(42)
for (m in 1:n_members) {
  pct <- round((m / n_members) * 100)
  cat(sprintf("\r   -> Generating member %2d/%2d [%-20s] %3d%%", 
              m, n_members, paste(rep("=", floor(pct / 5)), collapse = ""), pct))
  flush.console()
  
  # 1. Perturbazione delle Condizioni Iniziali (maggiore variabilità)
  X_test_noisy <- X_attuale + rnorm(length(X_attuale), mean = 0, sd = 0.12)
  X_test_scaled_noisy <- matrix((X_test_noisy - mean_X) / sd_X, nrow = 1)
  
  # 2. Diversificazione dell'Architettura della Rete (Size & Decay variabili)
  size_m  <- sample(8:16, 1)
  decay_m <- runif(1, 0.01, 0.2)
  
  set.seed(m * 105)
  mod <- nnet(X_train_scaled, Y_train_scaled, size = size_m, linout = TRUE, 
              maxit = 3000, trace = FALSE, decay = decay_m)
  
  pred_base <- as.numeric(predict(mod, newdata = X_test_scaled_noisy) * sd_Y + mean_Y)
  
  # 3. Iniezione di Rumore Stocastico Cumulativo (Random Walk settimanale)
  # Cresce nel tempo per dare la tipica forma a "spaghetto" frastagliato
  stochastic_noise <- cumsum(rnorm(n_ahead, mean = 0, sd = 0.04))
  
  ensemble_matrix_pure[m, ] <- pred_base + stochastic_noise
}
cat("\n   -> Ensemble completato!\n")

forzante_settimanale <- 0.65 * exp(-((1:n_ahead - 13)^2) / 18)
ensemble_matrix_forced <- t(apply(ensemble_matrix_pure, 1, function(x) x + forzante_settimanale))


# ==========================================
# 4. PLOTTING: DISPLAY A SCHERMO E SALVATAGGIO (CON EL NIÑO PRECEDENTE)
# ==========================================
cat("[4/4] Rendering dei grafici a schermo e salvataggio file in corso...\n")

if (!dir.exists("plots")) dir.create("plots")

# Estensione della serie storica a 160 settimane (~3 anni) per includere il picco 2023/2024
n_weeks_hist <- 160
date_storiche <- tail(dati_raw$Week, n_weeks_hist)
date_forecast <- seq(max(date_storiche) + weeks(1), by = "1 week", length.out = n_ahead)
valori_storici <- tail(anom_34, n_weeks_hist)

fun_plot_ensemble <- function(ens_mat, titolo, file_name) {
  all_dates <- c(date_storiche, date_forecast)
  n_hist <- length(date_storiche)
  val_t0 <- tail(valori_storici, 1)
  x_proj <- n_hist:(n_hist + n_ahead)
  
  ens_mean <- apply(ens_mat, 2, mean)
  ens_max  <- apply(ens_mat, 2, max)
  ens_min  <- apply(ens_mat, 2, min)
  
  y_lims <- range(c(valori_storici, ens_mat)) + c(-0.3, 0.3)
  
  disegna_grafico <- function() {
    plot(1:n_hist, valori_storici, type = "l", lwd = 2.5, col = "black",
         xlim = c(1, length(all_dates)), ylim = y_lims,
         main = titolo, xlab = "", ylab = "Anomalia SST NINO3.4 (°C)", xaxt = "n", las = 1)
    
    # Marcatori asse x ogni 12 settimane (~3 mesi) per mantenere la legenda leggibile
    at_ticks <- seq(1, length(all_dates), by = 12)
    axis(1, at = at_ticks, labels = format(all_dates[at_ticks], "%b %Y"), las = 2, cex.axis = 0.8)
    
    # Traiettorie dei singoli membri
    for (m in 1:nrow(ens_mat)) {
      lines(x_proj, c(val_t0, ens_mat[m, ]), col = rgb(0.15, 0.35, 0.75, 0.30), lwd = 1.1)
    }
    
    # Inviluppo Estremi e Media
    lines(x_proj, c(val_t0, ens_max), col = "firebrick", lty = 2, lwd = 2)
    lines(x_proj, c(val_t0, ens_min), col = "darkblue", lty = 2, lwd = 2)
    lines(x_proj, c(val_t0, ens_mean), col = "black", lwd = 3.5)
    
    # Linee di riferimento
    abline(h = 0, col = "darkgray", lty = 3)
    abline(h = 0.5, col = "red", lty = 2)
    abline(h = 1.5, col = "darkred", lty = 2)
    abline(v = n_hist, col = "gray40", lty = 3)
    
    legend("topleft", 
           legend = c("Dato Osservato", "Media Ensemble", "Scenari Estremi (Min/Max)", "Membri Ensemble (50)"),
           col = c("black", "black", "firebrick", rgb(0.15, 0.35, 0.75, 0.6)), 
           lwd = c(2.5, 3.5, 2, 1.1), lty = c(1, 1, 2, 1), bty = "n", cex = 0.85)
  }

  png(filename = file.path("plots", file_name), width = 2400, height = 1800, res = 300)
  disegna_grafico()
  dev.off()
  
  if (exists("dev.new")) dev.new()
  disegna_grafico()
}

fun_plot_ensemble(ensemble_matrix_pure, 
                  "18-Week Pure Dynamic Ensemble (50 Members, No Forcing)", 
                  "ensemble_pure_dynamic.png")

fun_plot_ensemble(ensemble_matrix_forced, 
                  "18-Week Hybrid Dynamic Ensemble (50 Members + Forcing)", 
                  "ensemble_forced_dynamic.png")

cat("   -> Elaborazione completata con successo!\n")
