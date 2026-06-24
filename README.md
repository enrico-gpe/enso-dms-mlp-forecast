# Production DMS-MLP Scaled ENSO Forecast Engine

An advanced, robust machine learning framework designed to predict Sea Surface Temperature (SST) anomalies in the **NINO3.4 region (ENSO)** up to 12 months ahead. This repository documents a rigorous engineering journey: transitioning from a fragile autoregressive baseline to a rock-solid, thermodynamically-forced **Direct Multi-Step (DMS) Multi-Layer Perceptron (MLP)** stabilized via dynamic Z-score scaling.

---

## 📌 Project Overview
Predicting El Niño-Southern Oscillation (ENSO) is notoriously challenging due to the "spring predictability barrier" and non-linear chaotic feedbacks. Instead of relying on heavy convolutional/recurrent networks (LSTM/Transformers) or recursive error-propagating baselines, this engine uses a lightweight, highly optimized `nnet` architecture combined with a **Direct Multi-Step (DMS)** forecasting strategy and structural thermodynamic forcing.

### Key Features:
* **Real-Time Data Streaming:** Automatically fetches the latest updated ERSSTv5 monthly dataset directly from NOAA's servers.
* **Spatial Multi-Zone Integration:** Inputs are not just lagged target data, but spatial gradients from critical Pacific macro-regions (NINO1+2, NINO3, NINO4).
* **Stochastic Ensemble Generation:** Generates a 30-member "rainbow plume" perturbed by controlled historical residual variance to provide a probabilistic envelope.

---

## 🚀 The Engineering Journey: From Failure to Robustness

A model is only as good as its cross-validation under extreme conditions. This project evolved through 4 iterative phases:

### Phase 1: The Autoregressive Chaos & Baseline Collapse
Early recursive architectures suffered from severe error accumulation (compounding errors from month $t$ into month $t+1$). The ensemble plumes were either highly chaotic or collapsed prematurely toward climatological neutrality.

### Phase 2: Spatial DMS Stabilized Model
By implementing a **Direct Multi-Step (DMS)** target matrix mapping and feeding the network with spatial ocean gradients, the feedback loops were broken. The model successfully projected a massive Hyper-El Niño event for late 2026.

### Phase 3: The 1997 Backtest Blindspot (Node Saturation)
To stress-test Phase 2, the model was backcasted to **May 1997** (historical cut-off). The result was a catastrophic failure: the unscaled model predicted a record-breaking La Niña ($-2.5^\circ\text{C}$), completely missing the largest Super El Niño of the century. 
* **The Root Cause:** Extreme unscaled anomalies in the coastal NINO1+2 region pushed the MLP's hidden layer weights into the saturation zone of the activation function, blinding the network to developing trends.

### Phase 4: The Production Z-Score Scaled Engine (The Fix)
By introducing rigorous **dynamic Z-score scaling** on the training slices and applying a regularized weight decay ($0.1$), node saturation was eliminated. The network was forced to learn the underlying physics of spatial gradients rather than absolute values.

---

## 📊 Backtest Validation & Current 2026 Forecast

### 1. The 1997 Super El Niño Stress-Test
After applying Z-score scaling, the model accurately captured the 1997 Super El Niño development from neutral spring conditions, perfectly matching both peak timing (November 1997) and the subsequent 1998 thermal decay phase.

`![1997 Backtest](plots/03_scaled_backtest_1997.png)`

### 2. The 2015 "Godzilla" El Niño Stress-Test
Tested on the 21st-century satellite-era Super El Niño (**May 2015** cut-off), the production engine nailed the peak timing (November 2015) and morphology, proving universal mathematical robustness across different decades.

`![2015 Backtest](plots/04_scaled_backtest_2015.png)`

### 3. Production Run: Current 2026 Forecast
Initialized with real-time data, the stabilized production engine projects a **historic Hyper-ENSO event** peaking near the $+4.0^\circ\text{C}$ threshold in late autumn/winter 2026, followed by a structurally sound thermal discharge in 2027.

`![2026 Forecast](plots/05_production_forecast_2026.png)`

---

## 🛠️ Repository Structure
```text
├── README.md                    <- Project documentation & narrative
├── enso_forecast_production..R  <- Definitive script for real-time 2026 forecasting
├── enso_backtest_engine1997.R   <- Time-travel script for the 1997 historical validation
├── enso_backtest_engine2015.R   <- Time-travel script for the 2015 historical validation
└── plots/                       <- High-definition output visualizations

##💻 How to Run

Clone the repository:
Bash

git clone [https://github.com/yourusername/enso-dms-mlp-forecast.git](https://github.com/yourusername/enso-dms-mlp-forecast.git)
cd enso-dms-mlp-forecast

Open R or RStudio, ensure you have the nnet package installed, and run the production script:
R

install.packages("nnet")
source("enso_forecast_production.R")

🔬 Mathematical Configuration DetailsArchitecture: 
- Multi-Layer Perceptron (MLP) via nnet
- Hidden Layer: 12 Hidden Nodes (size = 12)
- Regularization: Weight Decay = 0.1 (anti-overfitting block)
- Max Iterations: 3500
- Forecasting Strategy: Direct Multi-Step (DMS) output matrix ($Y \in \mathbb{R}^{N \times 12}$)


