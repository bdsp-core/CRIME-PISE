# CRIME_PISE

A machine learning model for predicting Post-Ischemic Stroke Epilepsy (PISE) and death in stroke patients using Random Survival Forest with competing risk analysis.

## Overview

CRIME_PISE is a prognostic tool that predicts the risk of unprovoked seizures (PISE) and mortality following ischemic stroke. The model uses clinical, imaging, and EEG features to generate individual risk scores and cumulative incidence functions (CIF) for competing outcomes.

The model was trained on 230 patients admitted between 2014-2020 and validated on 50 patients from 2021-2022.

## Features

The model incorporates the following features:

### Clinical Features
- **Age** at stroke onset
- **Pre-morbid status** (modified Rankin Scale)
- **Atrial fibrillation** history
- **NIHSS scores** (admission and 3-day post-stroke)

### SeLECT Score Components
- Stroke severity (NIHSS-based)
- Large artery atherosclerosis etiology
- Early seizures (≤7 days post-stroke)
- Cortical involvement
- MCA territory involvement

### Imaging Features
- Manual infarct volume measurement
- ACA territory involvement
- Cardioembolism etiology

### EEG Features
- **Epileptiform abnormality burden** (calculated using SPaRCNet)
- **Global theta power** and rhythmic activities
- **Hemispheric asymmetry** in power and rhythmic activities

## Repository Structure

```
CRIME_PISE/
├── code/
│   ├── CRIME_PISE.RData          # Pre-trained Random Survival Forest model
│   ├── example_code.qmd          # Example analysis workflow (Quarto)
│   ├── example_code.pdf          # Rendered example documentation
│   ├── funcs.R                   # Helper functions
│   ├── references.bib            # Bibliography
│   ├── data/
│   │   ├── data_train.csv        # Training dataset (230 patients)
│   │   └── data_test.csv         # Testing dataset (50 patients)
│   └── result/
│       ├── data_test_pred.csv    # Predictions on test set
│       └── data_test_eval.csv    # Evaluation metrics
└── README.md
```

## Installation

### Requirements

- R (≥ 4.0.0)
- Required R packages:
  - `dplyr`
  - `survival`
  - `prodlim`
  - `ggplot2`
  - `randomForestSRC`
  - `timeROC`
  - `shapper` (for SHAP analysis)

### Installing Dependencies

```r
install.packages(c("dplyr", "survival", "prodlim", "ggplot2",
                   "randomForestSRC", "timeROC"))
```

For `shapper` installation, follow instructions at: https://github.com/ModelOriented/shapper

## Usage

### 1. Load the Model and Data

```r
# Load helper functions
source("code/funcs.R")

# Load data
df.train <- read.csv('code/data/data_train.csv')
df.test <- read.csv('code/data/data_test.csv')

# Load pre-trained model
load("code/CRIME_PISE.RData")
```

### 2. Generate Risk Predictions

```r
# Predict PISE and Death risk scores
df.test$PISE_score <- predict(CRIME_PISE, df.test)$predicted[,1]
df.test$Death_score <- predict(CRIME_PISE, df.test)$predicted[,2]

# Stratify patients into risk groups
t.PISE <- quantile(df.train$PISE_score, 0.5)
t.Death <- quantile(df.train$Death_score, 0.5)

df.test$pred_label <- crime.predict(
  t.PISE = t.PISE,
  t.Death = t.Death,
  PISE.score = df.test$PISE_score,
  Death.score = df.test$Death_score
)
```

### 3. Individual Patient Analysis

```r
# Analyze first patient
i <- 1
obs_i <- df.test[i, CRIME_PISE$xvar.names]
pred_i <- predict(CRIME_PISE, obs_i)

# Plot cumulative incidence functions
plot_individual_cif(pred_i, df.test$event[i], df.test$time[i],
                    df.train$time, df.train$event)

# Generate SHAP values for feature importance
shap.df <- make_shap_df(mdl=CRIME_PISE, obs_i=obs_i)
plot_individual_shap(shap.df)
```

### 4. Model Evaluation

```r
# Evaluate model performance at 1-year post-stroke
eval.metrics <- crime.evaluate(
  yr = 1,
  t = df.test$time,
  y = df.test$event,
  yhat = df.test$pred_label,
  PISE.score = df.test$PISE_score,
  Death.score = df.test$Death_score
)
```

## Risk Stratification

Patients are stratified into four groups based on median risk scores:

1. **Event-Free**: Below-median for both PISE and Death
2. **PISE**: Above-median for PISE, below-median for Death
3. **Death**: Above-median for Death, below-median for PISE
4. **PISE or Death**: Above-median for both (further classified by higher score)

## Evaluation Metrics

The model reports:
- **AUC** for PISE and Death outcomes (time-dependent with competing risk)
- **Sensitivity** and **Specificity**
- **Positive/Negative Predictive Values** (PPV/NPV)
- **Event counts** among predicted positive cases

## Data Format

### Required Columns

- **Outcomes**: `event` (PISE/Death/Censored), `event_sd` (0/1/2), `time` (years)
- **Clinical**: `age`, `pre_morbid`, `afib`, `cardioembolism`, `nih_3d`
- **SeLECT Score**: `severity`, `large_artery`, `early_seizure`, `cortical_involvement`, `mca_involvement`, `select_score`
- **Imaging**: `hand_volume`, `aca`
- **EEG**: `total_ea_95`, `fft_theta_all_global`, `rhm_theta_all_global`, `fft_total_all_asym`, `rhm_total_all_asym`

## References

- **SeLECT Score**: Galovic et al. (2018). *The Lancet Neurology*, 17(2), 143-152. [DOI: 10.1016/s1474-4422(17)30404-0](http://dx.doi.org/10.1016/S1474-4422(17)30404-0)

- **SPaRCNet**: Ge et al. (2021). *Journal of Neuroscience Methods*, 351, 108966. [DOI: 10.1016/j.jneumeth.2020.108966](http://dx.doi.org/10.1016/j.jneumeth.2020.108966)

## Example Workflow

See [example_code.qmd](code/example_code.qmd) or [example_code.pdf](code/example_code.pdf) for a complete walkthrough.

## License

Please refer to the repository license file for usage terms.

## Contact

For questions or issues, please open an issue on this repository.
