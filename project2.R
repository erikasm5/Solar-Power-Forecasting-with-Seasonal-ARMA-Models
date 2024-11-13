
library(tidyverse)
library(lubridate)
library(ggplot2)

## 2 Predicting monthly solar power

solar_data <- read.csv("datasolar.csv")
solar_data$date <- make_date(solar_data$year, solar_data$month, 1)

#plotting time series:

ggplot(solar_data, aes(x = date, y = power)) +
  geom_line() +
  labs(title = "Solar Power Over Time",
       x = "Date",
       y = "Power (MWh)")


#####################################################################
#####################################################################
#####################################################################



#2.1

phi_1 <- -0.38
Phi_1 <- -0.94
mu <- 5.72
sigma_epsilon <- sqrt(0.222)  

solar_data <- solar_data %>%
  mutate(Xt = log(power) - mu)



### Vectorized solution
# Calculate Xhat
solar_data$Xhat2 <- mu + phi_1 * lag(solar_data$Xt, 1) + Phi_1 * lag(solar_data$Xt, 12)

# Calculate residuals
solar_data$residuals2 <- solar_data$Xt - solar_data$Xhat2
###



### For loop solution
Xhat <- numeric(nrow(solar_data))
residuals <- numeric(nrow(solar_data))

# Loop over time periods
for (i in 1:nrow(solar_data)) {
  if (i == 1) {
    # For the first observation, use the mean mu as the prediction
    Xhat[i] <- mu
  } else if (i <= 12) {
    # For the first 12 observations, use the mean mu plus the AR component
    Xhat[i] <- mu + phi_1 * (solar_data$Xt[i-1])
  } else {
    # For subsequent observations, use the mean mu plus both AR and seasonal AR components
    Xhat[i] <- mu + phi_1 * (solar_data$Xt[i-1]) + Phi_1 * (solar_data$Xt[i-12])
  }
  
  # Calculate residuals
  residuals[i] <- solar_data$Xt[i] - Xhat[i]
}

# Add Xhat and residuals to solar_data
solar_data$Xhat <- Xhat
solar_data$residuals <- residuals
###



# Do a model validation by checking the assumptions f i.i.d. errors
ggplot(solar_data, aes(x = date, y = residuals)) +
  geom_line() +
  labs(title = "Residuals Over Time",
       x = "Date",
       y = "Residuals")

# Histogram of Residuals
ggplot(solar_data, aes(x = residuals)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Residuals",
       x = "Residuals",
       y = "Frequency")

# Autocorrelation Function (ACF) Plot
#remove na values from xhat
solar_data2 <- na.omit(solar_data)
acf(solar_data$residuals, main = "Autocorrelation Function (ACF) of Residuals")

# Partial Autocorrelation Function (PACF) Plot
pacf(solar_data$residuals, main = "Partial Autocorrelation Function (PACF) of Residuals")

# Ljung-Box Test
ljung_box_test <- Box.test(solar_data$residuals, lag = 12, type = "Ljung-Box")
ljung_box_test







#####################################################################
#####################################################################
#####################################################################



#2.2


# add 12 following rows too the solar_data dataframe
solar_data2 <- solar_data %>%
  bind_rows(data.frame(date = seq(as.Date("2011-01-01"), by = "month", length.out = 12),
                       power = rep(NA, 12),
                       Xt = rep(NA, 12),
                       Xhat2 = rep(NA, 12),
                       residuals2 = rep(NA, 12)))




### Do we need to only use t=36, meaning we don't use seasonal part? or for seasonal part we also use t=36
#Manual model
solar_data2$Xhat[37] <- mu + phi_1 * solar_data2$Xt[36] + Phi_1 * solar_data2$Xt[25]

for (i in 38:nrow(solar_data2)) {
  solar_data2$Xhat[i] <- mu + phi_1 * solar_data2$Xt[i-1] + Phi_1 * solar_data2$Xt[i-12]
}
###


###ARIMA
arima_model <- arima(solar_data2$Xt[1:36], order = c(1, 0, 1), seasonal = list(order = c(0, 1, 1), period = 12))
arima_forecast <- predict(arima_model, n.ahead = 12)
arima_forecast$pred
solar_data2$Xhat[37:48] <- arima_forecast$pred
###



#creating a table of Predicted Xhat also include collumn for converted power values
table <- data.frame(date = solar_data2$date[37:48],
                    power = exp(solar_data2$Xhat[37:48] +mu ),
                    Xhat = solar_data2$Xhat[37:48])


#add power values to solar_data2
solar_data2$powerPred[37:48] <- exp(solar_data2$Xhat[37:48]+mu)
solar_data2$powerPred <- exp(solar_data2$Xhat + mu)

#plotting the predicted values with the original values in deiffernet color
ggplot(solar_data2, aes(x = date)) +
  geom_line(aes(y = power, color = "Original")) +
  geom_line(aes(y = solar_data2$powerPred, color = "Predicted")) +
  labs(title = "Solar Power Over Time",
       x = "Date",
       y = "Power (MWh)")

#####################################################################
#####################################################################
#####################################################################





#2.3





##### 3
library(astsa)

#3.1
# Simulate a (1, 0, 0) (0, 0, 0)12 ARIMA model with phi = 0.6
sim_data <- arima.sim(n=100, model = list(ar = c(0,0,0,0,0,0,0,0,0,0,0,0.7), ma= -0.4))
acf(sim_data, main="ACF of (1, 0, 0) × (0, 0, 0)12 Model")

#Plot
#par(mfrow = c(1, 1)) # Set plotting area to 1 row, 2 columns
acf(sim_data, main="ACF of (1, 0, 0) × (0, 0, 0)12 Model")
#pacf(sim_data, main="PACF of (1, 0, 0) × (0, 0, 0)12 Model")





#3.2
set.seed(123)
# Simulate a (0, 0, 0) (1, 0, 0)12 ARIMA model with seasonal_phi = -0.9
frequency <- 12 # Indicating monthly data with an annual cycle


model_spec <- Arima(ts(rnorm(100), frequency=frequency),
                    order=c(0,0,0), 
                    seasonal=list(order=c(1,0,0), period=frequency),
                    fixed=c(seasonal_phi=-0.9),
                    include.constant=FALSE)

# Simulate data from the model
sim_data <- simulate(model_spec, nsim=100)

# Plot the simulated data
plot(sim_data, main="Simulated SARIMA (1, 0, 0)x(1, 0, 0)[12] Data")
# Plot the ACF and PACF of the simulated data
par(mfrow=c(1, 2))
acf(sim_data, main="ACF of Simulated Data")
pacf(sim_data, main="PACF of Simulated Data")






#3.3
set.seed(123)
# Simulate a (1, 0, 0) (0, 0, 1)12 ARIMA model with phi = 0.9 and seasonal theta = -0.7
frequency <- 12 # Indicating monthly data with an annual cycle


model_spec <- Arima(ts(rnorm(100), frequency=frequency),
                    order=c(1,0,0), 
                    seasonal=list(order=c(0,0,1), period=frequency),
                    fixed=c(phi=-0.9, seasonal_theta=-0.7),
                    include.constant=FALSE)

# Simulate data from the model
sim_data <- simulate(model_spec, nsim=100)

# Plot the simulated data
plot(sim_data, main="Simulated SARIMA (1, 0, 0)x(1, 0, 0)[12] Data")
# Plot the ACF and PACF of the simulated data
par(mfrow=c(1, 2))
acf(sim_data, main="ACF of Simulated Data")
pacf(sim_data, main="PACF of Simulated Data")





#3.4
set.seed(123) # For reproducibility
# Simulate a (1, 0, 0) x (1, 0, 0)12 model with phi1 = -0.6 and Phi1 = -0.8
frequency <- 12 # Indicating monthly data with an annual cycle


model_spec <- Arima(ts(rnorm(100), frequency=frequency),
                    order=c(1,0,0), 
                    seasonal=list(order=c(1,0,0), period=frequency),
                    fixed=c(phi=-0.6, seasonal_phi=-0.8),
                    include.constant=FALSE)

# Simulate data from the model
sim_data <- simulate(model_spec, nsim=100)

# Plot the simulated data
plot(sim_data, main="Simulated SARIMA (1, 0, 0)x(1, 0, 0)[12] Data")
# Plot the ACF and PACF of the simulated data
par(mfrow=c(1, 2))
acf(sim_data, main="ACF of Simulated Data")
pacf(sim_data, main="PACF of Simulated Data")




#3.5
set.seed(123) # For reproducibility
# Simulate a (0, 0, 1) x (0, 0, 1)12 model with theta = 0.4 and seasonal_theta = -0.8
frequency <- 12 # Indicating monthly data with an annual cycle


model_spec <- Arima(ts(rnorm(100), frequency=frequency),
                    order=c(0,0,1), 
                    seasonal=list(order=c(0,0,1), period=frequency),
                    fixed=c(theta=0.4, seasonal_theta=-0.8),
                    include.constant=FALSE)

# Simulate data from the model
sim_data <- simulate(model_spec, nsim=100)

# Plot the simulated data
plot(sim_data, main="Simulated SARIMA (1, 0, 0)x(1, 0, 0)[12] Data")
# Plot the ACF and PACF of the simulated data
par(mfrow=c(1, 2))
acf(sim_data, main="ACF of Simulated Data")
pacf(sim_data, main="PACF of Simulated Data")




#3.6 
set.seed(123) # For reproducibility
# Simulate a (0, 0, 1) x (1, 0, 0)12 model with theta = -0.4 and seasonal_phi = 0.7
frequency <- 12 # Indicating monthly data with an annual cycle


model_spec <- Arima(ts(rnorm(100), frequency=frequency),
                    order=c(0,0,1), 
                    seasonal=list(order=c(1,0,0), period=frequency),
                    fixed=c(theta=-0.4, seasonal_theta=-0.7),
                    include.constant=FALSE)

# Simulate data from the model
sim_data <- simulate(model_spec, nsim=100)

# Plot the simulated data
plot(sim_data, main="Simulated SARIMA (1, 0, 0)x(1, 0, 0)[12] Data")
# Plot the ACF and PACF of the simulated data
par(mfrow=c(1, 2))
acf(sim_data, main="ACF of Simulated Data")
pacf(sim_data, main="PACF of Simulated Data")






ts_data <- ts(solar_data$power, frequency=12)
model_spec <- Arima(ts_data,
                    order=c(0,0,1), 
                    seasonal=list(order=c(1,0,0), period=frequency),
                    fixed=c(theta=0.5, seasonal_theta=-0.7),
                    include.constant=FALSE)

# Simulate data from the model
sim_data <- simulate(model_spec, nsim=100)

# Plot the simulated data
plot(sim_data, main="Simulated SARIMA (1, 0, 0)x(1, 0, 0)[12] Data")
# Plot the ACF and PACF of the simulated data
par(mfrow=c(1, 2))
acf(sim_data, main="ACF of Simulated Data")
pacf(sim_data, main="PACF of Simulated Data")









#Recalculted


#3.6
# Parameters
theta <- -0.4
Phi <- 0.7

# Autocorrelation calculations
rho <- numeric(36) # Calculate ACF for 36 lags as an example
rho[1] <- 1 # rho(0)
rho[2] <- theta / (1 + theta^2) # rho(1)
rho[12] <- theta * Phi / (1 + theta^2) # rho(11), since R indexing starts at 1

# Calculate ACF for k >= 12 using recursive formula
for (k in 13:36) {
  rho[k] <- Phi * rho[k-12]
}

# Lags
lags <- 0:35

# Plot ACF
plot(lags, rho, type="h", lwd=2, col="blue", xlab="Lag", ylab="Autocorrelation", main="ACF for (0, 0, 1)x(1, 0, 0)[12] Model")
abline(h=0, col="red", lty=2)






set.seed(123) # For reproducibility

# Define the coefficients
ma_coeffs <- c(-0.4) # Non-seasonal MA(1) parameter
ar_coeffs <- rep(0, 12) # Initialize AR parameters for 12 lags
ar_coeffs[12] <- 0.7 # Set the seasonal AR(1) parameter at lag 12

# Simulate the time series
n <- 120 # Length of the simulated time series
sim_data <- arima.sim(n = n, model = list(order = c(0, 0, 1), ma = ma_coeffs, ar = ar_coeffs))

# Plot the simulated time series
plot(sim_data, main="Simulated Time Series: (0, 0, 1)x(1, 0, 0)[12]")

# Plot the ACF and PACF
par(mfrow=c(1, 2))
acf(sim_data, main="ACF of Simulated Time Series")
pacf(sim_data, main="PACF of Simulated Time Series")



library(forecast)
set.seed(111)
#3.1
set.seed(111)
# Simulate a (1, 0, 0) (0, 0, 0)12 ARIMA model with phi = 0.6
sim_data <- arima.sim(n=100, model = list(ar = c(-0.6)))
acf(sim_data, main="ACF of (1, 0, 0) × (0, 0, 0)12 Model",lag = 100)
pacf(sim_data, main="PACF of (1, 0, 0) × (0, 0, 0)12 Model",lag = 100)
plot(sim_data)


#3.2
set.seed(111)
# Simulate a (0, 0, 0) (1, 0, 0)12 ARIMA model with seasonal_phi = -0.9
sim_data <- arima.sim(n=100, model = list(ar = c(0,0,0,0,0,0,0,0,0,0,0,0.9)))
acf(sim_data, main="ACF of (0, 0, 0) (1, 0, 0)12 Model",lag = 100)
pacf(sim_data, main="PACF of (0, 0, 0) (1, 0, 0)12 Model",lag = 100)
plot(sim_data)



#3.3
set.seed(111)
# Simulate a (1, 0, 0) (0, 0, 1)12 ARIMA model with phi = 0.9 and seasonal theta = -0.7
sim_data <- arima.sim(n=100, model = list(ar = c(-0.9), ma= c(0,0,0,0,0,0,0,0,0,0,0,-0.7)))
acf(sim_data, main="ACF of (1, 0, 0) (0, 0, 1)12 Model",lag = 100)
pacf(sim_data, main="PACF of (1, 0, 0) (0, 0, 1)12 Model",lag = 100)
plot(sim_data)



#3.4
set.seed(111)
# Simulate a (1, 0, 0) x (1, 0, 0)12 model with phi1 = -0.6 and Phi1 = -0.8


# A simulation function for ARMA simulation, use model as arima.sim, i.e. flip sign of phi (into ar) coefficients
sim <- function(model, n=100, nburnin=100){
  n <- n + nburnin
  # Take the ar and ma part
  ar <- model$ar
  ma <- model$ma
  # The order (i.e. the number of lags)
  p <- length(ar)
  q <- length(ma)
  # The vector for the simulation result
  y <- numeric(n)
  # Generate the random normal values
  eps <- rnorm(n)
  # Run the simulation
  for(i in (max(p,q)+1):n){
    y[i] <- eps[i] + sum(y[i-(1:p)] * ar)
  }
  # Return without the burn-in period
  return(y[(nburnin+1):n])
}


# Non-stationary process
# Do the simulation and plot
n <- 100
model <- list(ar = c(0.6, 0,0,0,0,0,0,0,0,0,0, 0.8, 0.8*0.6))
sim_data <- sim(model, n)


# Plotting the ACF of the simulated series
acf(sim_data, main="ACF of Simulated (1, 0, 0)x(1, 0, 0)[12] Model",lag = 100)
pacf(sim_data, main="PACF of Simulated (1, 0, 0)x(1, 0, 0)[12] Model",lag = 100)
plot(sim_data, type="l", ylab="x", xlab="t")




library(sarima)
library(portes)
# Question 3
# (mixed) seasonal model (p, 0, q)×(P, 0, Q)
set.seed(111)
n <- 100

x <- sim_sarima(n=n,model=list(ar=0.6,sar =0.8, nseasons=12, sigma2 = 1))

acf(x,lag.max = 100, main=NA)  
pacf(x,lag.max = 100, main=NA)
plot(x,ylab= expression(Y[t]), xlab = "Time", type = "l")


sim <- function(n, phi, Phi, nburnin=100) {
  n_total <- n + nburnin  # Adjust for burn-in period
  y <- rep(0, n_total)    # Initialize the series with zeros
  
  # Pre-generate the white noise (epsilon terms)
  eps <- rnorm(n_total)
  
  # Simulation must consider the burn-in period to ensure stability
  for(i in 2:n_total) {
    # Directly apply non-seasonal AR component for all
    temp <- phi * y[i - 1] + eps[i]
    
    # Apply seasonal AR component if beyond the first seasonal period
    if(i > 12) {
      temp <- temp + Phi * y[i - 12]
    }
    
    # Assign the calculated value back to the series
    y[i] <- temp
  }
  
  # Return the series excluding the burn-in period
  return(y[(nburnin + 1):n_total])
}






#3.5
set.seed(111)
# Simulate a (0, 0, 1) x (0, 0, 1)12 model with theta = 0.4 and seasonal_theta = -0.8
sim_data <- arima.sim(n=100, model = list(ma= c(0.4, 0,0,0,0,0,0,0,0,0,0, -0.8, -0.8*0.4)))
acf(sim_data, main="ACF of (0, 0, 1) x (0, 0, 1)12 Model",lag = 100)
pacf(sim_data, main="PACF of (0, 0, 1) x (0, 0, 1)12 Model",lag = 100)
plot(sim_data)


#3.6
set.seed(100)
# Simulate a (0, 0, 1) x (1, 0, 0)12 model with theta = -0.4 and seasonal_phi = 0.7
sim_data <- arima.sim(n=100, model = list(ar = c(0,0,0,0,0,0,0,0,0,0,0,-0.7), ma= -0.4))
acf(sim_data, main="ACF of (0, 0, 1) x (1, 0, 0)12 Model", lag = 100)
pacf(sim_data, main="PACF of (0, 0, 1) x (1, 0, 0)12 Model",lag = 100)
plot(sim_data)


