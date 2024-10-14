########################################################################
#### Introduction ####
# We will analyze the problem of Long-Run investing in the case where the investment opportunity set (IOS) is constant,
# and there is no intermediate consumption/further investment into the portfolio.
# We will assume CRAA risk-preference and we want to make distributional assumsption.
# We will assume that returns are jointly log-normally distributed, this is a more appropiate compared to normality.
# Reason for this is the following:
  # 1. Log-normality does not allow for gross returns falling below -100%
  # 2. Compounded returns are also log-normality distributed because the log-returns add up over time


########################################################################
#### Preparatory Steps ####

## To clear the working station
rm(list = ls())

username = Sys.getenv("USERNAME")

## Setting the working directory
setwd(paste0("C:\\Users\\", username, "\\OneDrive - BI Norwegian Business School (BIEDU)\\Desktop\\BI-Skole\\3 Semester\\Strategic Asset\\Homework_tasks"))
source(".\\2\\portConstruct.R")
source(".\\3\\portUtility.R")
source(".\\3\\variousTools.R")
library("MASS")
library("moments")


######################### Import data and explanation #########################
data = read.csv(".\\1\\test_data.csv")

# Choose risky asset set
asset_set = 1
if(asset_set == 1){
  # We choose bonds and stocks
  data_ret = data[, c(2, 4)] / 100;
} else if (asset_set == 0){
  # Bonds and 5 industries
  data_ret = data[, 4:9] / 100;
}


# Excess returns
data_rf = data[, 3] / 100
data_eret = sweep(data_ret, 1, data_rf, "-")

# Mean and variance of one-period excess returns
Mu_e = colMeans(data_eret)
Sigma = cov(data_eret)

# Number of risky assets
n_assets = length(Mu_e)


######################### 1a. Optimal static portfolios #########################
# We can try the following approaches for estimating the optimal portfolio
# 1. Mean-variance investing with one period parameter estimates
# 2. Mean-variance investing with T-horizon parameter estimates
# 3. Expected utility maximization using log-normality assumption
# 4. Expected utility maximization using historical simulation
# 3 and 4 are the correct approaches and 1 and 2 are approximations. So let us start with computing optimal portfolios with 1 and 2. 


# Horizon (10 year = 120 months)
T = 120

# Relative risk-aversion
gamma = 5

# Risk free rate
yT = 0.001
Rf_T = (1 + yT)^T

# Mean and variance of gross returns of risky assets
Mu = 1 + yT + colMeans(data_eret)
Sigma = cov(data_eret)

# T horizon meand and variance of gross returns
Mu_T = Mu^T

# COmpute covariance matrix for horizon T with iterations
Mu2 = as.matrix(Mu) %*% Mu
Sigma_T = Sigma

for(t in seq(2,T)){
  Sigma_T = Sigma_T * (Sigma + Mu2) + Sigma * (Mu2^(t-1))
}

# One period optimal portfolios without (w1) and with (w1f) risk-free asset
w1 = uncMeanVar_ra(Mu, Sigma, gamma)
w1t = uncTangent(Mu_e, Sigma)   # Tangency portfolio
l = c((w1t %*% Mu_e) / (gamma * (w1t %*% Sigma %*% w1t)))   # Optimal leverage
w1f = l * w1t

# T period optimal portfolios
wT = uncMeanVar_ra(Mu_T, Sigma_T, gamma)
wTt = uncTangent(Mu_T - Rf_T, Sigma_T)
l = c((wTt %*% (Mu_T - Rf_T)) / (gamma * (wTt %*% Sigma_T %*% wTt)))   # Optimal leverage
wTf = l * wTt

# Output results
t(rbind(w1, w1f))
t(rbind(wT, wTf))


## We will now solve the expected utility maximization problem assuming:
# 1. Log-normal distributed returns
# 2. No parametric assumption and historical returns are IID

# Log normality assumption without risk-free asset
optPort = portCRRA_LT(data = 1 + yT + data_eret, model = "logn", T = T, gamma = gamma,
                      N = 15000)
wTl = optPort$w

# Log normality assumption with risk-free asset 
optPort = portCRRA_LT(data = 1 + yT + data_eret, model = "logn", T = T,rf = log(1+yT) ,gamma = gamma,
                      minlev = -1, maxlev = 2, lower = 0, upper = 2, N = 15000)

wTlf = optPort$w

# Historical simulations without risk-free asset
optPort <- portCRRA_LT(data=1 + yT + data_eret, model="hist", T=T, gamma=gamma, N=10000)
wTh <- optPort$w

# Historical simulations with risk-free asset
optPort <- portCRRA_LT(data=1 + yT + data_eret, model="hist", T=T, rf=log(1+yT), gamma=gamma, minlev=-1,
                       maxlev=2, lower=0, upper=2, N=10000)
wThf <- optPort$w

#output results
t(rbind(wTl,wTh))
t(rbind(wTlf,wThf))


######################### 1b. Optimal portfolios and investment horizon #########################
# We will vary the horizon, compute optimal static portfolios, assuming constant IOS (returns are IID)


# Various horizons in months:
T = c(1, 3, 6, 12, 24, 60, 120, 180, 240)
nT = length(T)

# We create arrays to store portfolio weights and CE
WW = array(0, c(nT, n_assets))
CE = array(0, c(nT, 1))

for (i in 1:nT){
  # Perform optimization for each horizon T
  optPort = portCRRA_LT(data = 1 + yT + data_eret, model = "logn", T = T[i],rf = NA, gamma = 5,
                        N = 15000)
  
  # Collect result: portfolio weights and CE
  WW[i, ] = optPort$w
  CE[i] = optPort$ce 
  
  }

# Plot the weight of the stock (or/and other portfolio weights) against horizon
if (asset_set == 1){
  # Plot the stock weight
  plot(T, WW[,1], type="o", cex=1.5, lwd=2, col="blue", main="Stocks and bonds", xlab="Horizon (months)", 
       ylab="Stock weight", cex.lab=1.5, cex.main=2, cex.axis=1.2)
} else {
  # Plot all weights
  matplot(T, WW, type = "b",pch=1,cex=1.5,lwd=2,main="Bonds and 5 industries",ylab="Weights",
          xlab="Horizon (months)", cex.lab=1.5, cex.main=2, cex.axis=1.2)
  legend("topright", legend = c("bonds","ind1","ind2","ind3","ind4","ind5"), col=1:6, pch=1)
}

###################################################################################################

######################### 2. Rebalancing VS Buy-and-Hold #########################

# In this section we compare the rebalancing and the buy-and-hold strategies.
# We compare the optimal static long-run portfolio for the buy-and-hold strategy with the static one period optimal portfolio rebelanced over time
# We also use the Geometric Mean in order to compare total return of the strategy for a certain path

# Horizon
T = 120

# Number of simulation paths
N = 15000
set.seed(1)

# Choose BnH: BnH_short = TRUE means that it starts with the one-period optimal portfolio
# BnH_short = FALSE means that it is the long-term optimal portfolio
bnh_short = FALSE

# Optimal one period portfolio: rebalancing strategy
optPort = portCRRA_LT(data = 1 + yT + data_eret, model = "logn", T =1, rf=NA, gamma = gamma)
W_reb = as.matrix(optPort$w)

# Buy and hold strategy
if (bnh_short == TRUE){
  # Optimal short term portfolio
  W_bnh <- W_reb
} else {
  # Optimal static long-run portfolio
  optPort <- portCRRA_LT(data=1 + yT + data_eret, model="logn", T=T, rf=NA, gamma=gamma)
  W_bnh <- as.matrix(optPort$w)
}

# Arrays to store simulated portfolio total returns
r_reb <- rep(0, times=N)
r_bnh <- rep(0, times=N)

# Compute moments of log returns
Mu_ln = colMeans(log(1 + yT + data_eret))
Sigma_ln = cov(log(1 + yT + data_eret))

# Simulate N paths
for(i in 1:N){
  # Asset returns
  R = exp(mvrnorm(T, Mu_ln, Sigma_ln))
  r = log(R)
  
  # Rebelanced portfolio log total return
  r_reb[i] = sum(log(R %*% W_reb))
  
  # Buy-and-hold portfolio log total return
  r_bnh[i] = log(exp(colSums(r)) %*% W_bnh)
}

# Plot the two distributions of total returns
ecdf_reb <- ecdf(exp(r_reb))
ecdf_bnh <- ecdf(exp(r_bnh))
plot(ecdf_reb, verticals=TRUE, do.points=FALSE, main="Return distributions", col="blue", xlab="Total return", cex.lab=1.5, cex.main=2)
plot(ecdf_bnh, verticals=TRUE, do.points=FALSE, add=TRUE, col="red")
legend(x = "topleft", legend = c("Rebalanced", "Buy-and-hold"), lty = c(1, 1), col = c("blue", "red"), lwd = 1) 

# Plot the two distributions of geometric means
gm_reb <- 100*(exp(r_reb / T) - 1)
gm_bnh <- 100*(exp(r_bnh / T) - 1)
ecdf_g_reb <- ecdf(gm_reb)
ecdf_g_bnh <- ecdf(gm_bnh)
plot(ecdf_g_reb, verticals=TRUE, do.points=FALSE, main="Return distributions", col="blue", xlab="Geometric mean (%)", cex.lab=1.5, cex.main=2)
plot(ecdf_g_bnh, verticals=TRUE, do.points=FALSE, add=TRUE, col="red")
legend(x = "topleft", legend = c("Rebalanced", "Buy-and-hold"), lty = c(1, 1), col = c("blue", "red"), lwd = 1) 

# Statistics of the geometric means of the two strategies
Stats <- matrix(0, nrow=6, ncol=2)
colnames(Stats) <- c("Rebalanced", "Buy-and-hold")
rownames(Stats) <- c("Mean", "Median", "Std", "Skew", "Kurt", "CE")
Stats[1,] <- c(mean(gm_reb), mean(gm_bnh))
Stats[2,] <- c(median(gm_reb), median(gm_bnh))
Stats[3,] <- c(sd(gm_reb), sd(gm_bnh))
Stats[4,] <- c(skewness(gm_reb), skewness(gm_bnh))
Stats[5,] <- c(kurtosis(gm_reb), kurtosis(gm_bnh))

# Certainty equivalents
Stats[6,] <- c(mean(exp((1-gamma)*r_reb))^(1/(1-gamma)), mean(exp((1-gamma)*r_bnh))^(1/(1-gamma)))

options(digits=3)
Stats


# To compute the expectation we use numerical integration with the Gauss-Hermite quadrature points.
# This is a much more efficient method compared to Monte Carlo simulations
pts <- mgauss.hermite(20, mu=Mu_ln*T, sigma=Sigma_ln*T, prune=0.1)
R <- exp(pts$points)                 # Points of the distribution of returns
p <- pts$weights                     # Probability weights of the points

# Define the function to compute the CRRA certainty equivalent (CE)
certeq <- function(theta){
  util <- (R %*% theta)^(1-gamma)
  (p %*% util)^(1/(1-gamma))
}

# Choose the set of weights of the first asset to perform the CE computations
w <- seq(from=0.4, to=0.9, by=0.05)
n <- length(w)
CE <- rep(0, times=n)
for (i in 1:n){
  CE[i] <- certeq(c(w[i],1-w[i]))
}

# Plot CE against the weights
plot(w, CE, type="o", cex=1.5, lwd=2, col="blue", main=sprintf("Buy-and-hold (T=%2.0f)",T), xlab="Stock weight", ylab="Certainty equivalent", cex.lab=1.5,
     cex.main=2, cex.axis=1.2)
