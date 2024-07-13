# Assignment 3 sample R code

library(nloptr)
library(quantmod)
library(xts)
library(zoo)
library(data.table)
library(ggplot2)
library(MASS)
library(metRology)
library(moments)
library(parallel)
library(Quandl)
library(rugarch)
library(dplyr)

# F is the price of the underlying asset
# X is the strike price of the option
# t is the time to maturity (in years)
# r is the tbill rate (in decimal form)
# sigma is the volatility of the underlying asset
# BFC is the price of the call option
# BFP is the price of the put option
# IVC is the implied volatility of the call option
# IVP is the implied volatility of the put option

BFC <- function(F,X,t,r,sigma) {
  d1 <- log(F/X) + 0.5*sigma^2 *t
  d1 <- d1/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  N1 <- pnorm(d1)
  N2 <- pnorm(d2)
  C <- exp(-r*t) * (F * N1 - X * N2 )
  return(C)
}

BFP <- function(F,X,t,r,sigma) {
  d1 <- log(F/X) + 0.5*sigma^2 *t
  d1 <- d1/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  NM1 <- pnorm(-d1)
  NM2 <- pnorm(-d2)
  P <- exp(-r*t) * (X * NM2 - F * NM1 )
  return(P)
}

IVC <- function(F,X,t,r,Call) {
  eval_f_C <- function(sigma) {
    return ( (Call-BFC(F,X,t,r,sigma))^2 ) 
  }
  opts <- list("algorithm"="NLOPT_LN_COBYLA",
               "xtol_rel"=1.0e-8)
  xs <- 0.10
  es <- nloptr( x0=xs,
                eval_f=eval_f_C,
                opts=opts)
  return(es$solution)
}
IVP <- function(F,X,t,r,Put) {
  eval_f_P <- function(sigma) {
    return ( (Put-BFP(F,X,t,r,sigma))^2 ) 
  }
  opts <- list("algorithm"="NLOPT_LN_COBYLA",
               "xtol_rel"=1.0e-8)
  xs <- 0.10
  es <- nloptr( x0=xs,
                eval_f=eval_f_P,
                opts=opts)
  return(es$solution)
}

DATE <- function(yyyy,mm,dd) {
  dte  <- as.Date(sprintf("%i-%i-%i",yyyy,mm,dd),format="%Y-%m-%d")
  return(dte)
}

#----------------------------------------------------
setwd("/Users/zhangdonglei/Downloads")
load("MQM530-TeamAssignment1.RData", verbose=TRUE)

#############################################################
#Part 1
#----------------------------------------------------


load("MQM530-TeamAssignment1.RData")
# index(logret_X.xts)
# coredata(logret_X.xts)
logret = as.vector(coredata(logret_X.xts))
logret

# normal
fit_n <- fitdistr(logret,'normal')
mu <- fit_n$estimate[1]
sig<- fit_n$estimate[2]
AIC_n <- 2*length(fit_n$estimate) - 2*fit_n$loglik
cat("Mean: ",round(mu,6),'\n')
cat("SD: ",round(sig,6),'\n')
cat("Loglik:",round(fit_n$loglik,2),'\n')
cat("AIC_n:",round(AIC_n,2),'\n')

# scaled t
fit_t <- fitdistr(logret,'t')
round(fit_t$estimate,6)
round(fit_t$loglik,2)
m <- fit_t$estimate[1]
s <- fit_t$estimate[2]
tdf <- fit_t$estimate[3]

AIC_t <- 2*length(fit_t$estimate) - 2*fit_t$loglik
cat("AIC of scaled t = ",AIC_t,'\n')    
AIC_n > AIC_t

# Since the Akaike Information Criterion of a normal distribution model is greater than the one of a scaled-t distribution model, we should pick the latter model, scaled-t, for minimizing the AIC to get the best fit distribution. 

CalcVaRES <- function(r,alpha) {
  VaR <- quantile(r,1-alpha)
  ES  <- mean(r[r<VaR])
  VaR_ES <- c(VaR,ES)
  names(VaR_ES) <- c("VaR","ES")
  return(VaR_ES)
}

nsim <- 100000
nper <- 5
sim <- rep(0,nsim)
set.seed(437)
for (i in 1:nper) {
  sim <- sim+sample(logret,nsim,replace=TRUE)
}
alpha <- 0.95
VaR_ES <- CalcVaRES(sim,alpha)
round(VaR_ES,6)
      


#----------------------------------------------------

data <- copy(holdings_X)
setnames(data, "FutPrice" , "F")
setnames(data, "Strike" , "X")
data[, IVP := NA_real_] 
data[, IVC := NA_real_]


for (i in 1:nrow(data)) {
  row <- data[i]
  F <- row$F
  X <- row$X
  t <- as.numeric(difftime(row$Expiration, row$Date, units = "days")) / 365
  r <- row$r
  if (row[["Call/Put"]] == 1) {
    Call <- row$OptPrice
    data[i, IVC := IVC(F, X, t, r, Call)]
  } else if (row[["Call/Put"]] == 2) {
    Put <- row$OptPrice
    data[i, IVP := IVP(F, X, t, r, Put)]
  }
}

# checking patterns in strike price and implied volatility
data$IV <- coalesce(data$IVC, data$IVP)

ggplot(data, aes(x = X, y = IV)) +
  geom_line() +
  labs(x = "Strike Price", y = "Implied Volatility", title = "Overall Volatility Skewness")

ggplot(data, aes(x = X, y = IVP)) +
  geom_line() +
  labs(x = "Strike Price", y = "Implied Volatility", title = "Put Volatility Skewness")

ggplot(data, aes(x = X, y = IVC)) +
  geom_line() +
  labs(x = "Strike Price", y = "Implied Volatility", title = "Call Volatility Skewness") +
  xlim(2550, max(data$X))

#----------------------------------------------------

cash <- 628226078
# Function to calculate market value
calculate_market_value <- function(optprice, contracts) {
  return(optprice * contracts * 250)
}

# Given data
F <- data$F
X <- data$X
t <- as.numeric(data$Expiration - data$Date) / 365 
r <- data$r

# Calculate market value change for different IV multipliers
IV_multipliers <- c(1.5, 2, 3, 4)
contracts <- data$Contracts

for (i in 1:length(F)) {
  orig_IV <- data$IV[i]
  for (multiplier in IV_multipliers) {
    sigma <- multiplier * orig_IV
    optprice_call <- BFC(F[i], X[i], t[i], r[i], sigma)
    optprice_put <- BFP(F[i], X[i], t[i], r[i], sigma)
    market_value_call <- calculate_market_value(optprice_call, contracts[i])
    market_value_put <- calculate_market_value(optprice_put, contracts[i])
    market_value_change_call <- market_value_call / calculate_market_value(BFC(F[i], X[i], t[i], r[i], orig_IV), contracts[i])
    market_value_change_put <- market_value_put / calculate_market_value(BFP(F[i], X[i], t[i], r[i], orig_IV), contracts[i])
    cat("For row", i, "and IV multiplier:", multiplier, "\n")
    cat("Market value change percentage for Call option:", market_value_change_call, "\n")
    cat("Market value change percentage for Put option:", market_value_change_put, "\n\n")
  }
}

###########
# Given data
F <- data$F
X <- data$X
t <- as.numeric(data$Expiration - data$Date) / 365  # Convert to numeric and scale to years
r <- data$r
IV <- data$IV
contracts <- data$Contracts

# Calculate initial fund value
initial_value <- sum(calculate_value(BFC(F, X, t, r, IV), contracts),
                     calculate_value(BFP(F, X, t, r, IV), contracts))

# Define target value (50% reduction)
target_value <- initial_value * 0.5

# Incrementally increase IV until fund value decreases by 50%
new_IV <- IV
while (TRUE) {
  # Calculate new fund value
  new_value <- sum(calculate_value(BFC(F, X, t, r, new_IV), contracts),
                   calculate_value(BFP(F, X, t, r, new_IV), contracts))
  
  # Print values for debugging
  cat("IV:", new_IV, "\tNew value:", new_value, "\n")
  
  # Check if new fund value is below target value
  if (new_value <= target_value) {
    break
  }
  new_IV <- new_IV + 0.001  # Increment IV by 0.001 (adjust this increment as needed)
}

IV_increase <- new_IV - IV


#############################################################


uspec_t <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
                      distribution.model = "std")
ncore <- max(1, detectCores(logical = FALSE) / 2)
cl <- makeCluster(ncore)
n2020 <- length(logret_G.xts["/2020-12-31"])
logretg.xts <- logret_G.xts["/2021-06-30"]
roll_garch <- ugarchroll(spec = uspec_t, 
                         data = logret_G.xts, 
                         n.start = n2020, 
                         refit.every = 1, 
                         refit.window = "recursive", 
                         calculate.VaR = TRUE, 
                         VaR.alpha = 0.05, 
                         n.ahead = 1, 
                         cluster = cl,
                         solver.control = list(trace=0),
                         keep.coef = TRUE)
stopCluster(cl)
df <- roll_garch@forecast$VaR
names(df) <- c("VaR")
df$date <- index(logret_G.xts)[(n2020 + 1):length(logret_G.xts)]
position <- 1000000
MarginRequirement <- position*(exp(-df$VaR)-1)
df_margin <- data.table(Date = df$date, MarginRequirement = MarginRequirement)
print(df_margin)


logretG <- as.numeric(logret_G.xts)
spec <- ugarchspec(mean.model = list(armaOrder = c(0,0)),
                   variance.model = list(model = "sGARCH"),
                   distribution.model = "std")
VaR_data <- numeric(length(logretG))
margin <- numeric(length(logretG))
money <- 1000000
start <- as.Date("2021-01-01")
endD <- as.Date("2021-06-30")
for (i in seq_along(logretG)) {
  if (index(logret_G.xts)[i] >= start & index(logret_G.xts)[i] <= endD) {
    fit <- ugarchfit(spec, data = logretG[1:i])
    forecast <- ugarchforecast(fit, n.ahead = 1)
    VaR <- as.numeric(quantile(sigma(forecast), probs = 0.05))  
    VaR_data[i] <- VaR
    margin[i] <- money * (exp(VaR) - 1)
  }
}

actualret <- exp(logretG) - 1
below_VaR <- sum(actualret[index(logret_G.xts) >= start & index(logret_G.xts) <= endD] < VaR_data[index(logret_G.xts) >= start & index(logret_G.xts) <= endD])
print(below_VaR)
