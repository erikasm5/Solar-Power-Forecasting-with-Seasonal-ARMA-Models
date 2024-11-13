# Solar Power Forecasting with Seasonal ARMA Models

## Project Overview
This project focuses on predicting monthly solar power generation using seasonal ARMA models. The goal is to apply statistical time series methods to forecast energy production from a solar PV plant, leveraging historical data. This analysis helps in planning the operation of renewable energy systems by providing insight into both short-term and long-term power generation trends.

## Project Files
- **project2.R**: Main R script containing the analysis and forecasting code. The script:
  - Loads and preprocesses solar power data.
  - Constructs a seasonal ARMA model for monthly power generation.
  - Provides predictions for the next 12 months, along with prediction intervals.
  - Validates model assumptions and plots relevant time series and autocorrelation functions.

- **datasolar.csv**: Contains historical monthly solar power generation data used in the analysis.

- **assignment2_2024.pdf**: A guideline document describing the objectives and requirements for this analysis.

## Instructions
1. Ensure all required R packages (`tidyverse`, `lubridate`, `ggplot2`) are installed:
   ```r
   install.packages(c("tidyverse", "lubridate", "ggplot2"))

2. Follow assiignment2.pdf instructions and complementary project2.R code
   
