########################################################################
#### Preparatory Steps ####

## To clear the working station
rm(list = ls())

username = Sys.getenv("USERNAME")

## Setting the working directory
setwd(paste0("C:\\Users\\", username, "\\OneDrive - BI Norwegian Business School (BIEDU)\\Desktop\\BI-Skole\\3 Semester\\Strategic Asset\\Homework_tasks"))
source(".\\2\\portConstruct.R")

# The portconstruct.r contains the following functions:
# 1. portMinVar: Returns the unconstrained or constrained minimum variance portfolio
# 2. uncTangent: Returns the unconstrained mean-variance tangent portfolio
# 3. conTangent: Returns the constrained mean-variance tangent portfolio
# 4. uncMeanVar_ra: Unconstrained mean-variance optimal portfolio for given risk aversion parameter "gamma"
# 5. uncMeanVar_mu: Unconstrained mean-variance optimal portfolio for a target mean
# 6. uncMeanVar_sd: Unconstrained mean-variance optimal portfolio for a target standard deviation
# 7. conMeanVar: Constrained mean-variance optimal portfolio for a given risk-aversion, or for target mean, or for a target volatility
# 8. naiveVolRiskBudget: Returns portfolio weights based on naive volatility risk budgeting
# 9. naiveVarRiskBudget: Returns portfolio weights based on naive variance risk budgeting
# 10. trueVolRiskBudget: Returns portfolio weights based on true risk budgeting with respect to volatility contributions
# 11. maxDiversify: Returns portfolio weights of the unconstrained maximum diversification portfolio

######################### Import data and explanation #########################
# Load a test data set to perform some analysis. The data contains the time series of the following monthly returns in percentages:
# - Mkt: Market return
# - Rf: Risk-free return
# - Bnds10: Return on 10-year treasure bonds
# - Cnsmr: Return on consumer goods industry portfolio
# - Manuf: Return on manufacturing goods industry portfolio
# - HiTec: Return on high-tech industry portfolio
# - Hlth: Return on health industry portfolio
# - Other: Return on other industry portfolio
# - SMB: Return on SMB factor portfolio
# - HML: Return on HML factor portfolio
# - MOM: Return on momentum factor portfolio

# The first column contains the date of the returns.

data = read.csv(".\\1\\test_data.csv")


######################### Creation of portfolio #########################
# We will consider portfolios of the market (the five industries), and the bonds (with and WO a risk-free rate)
# We therefore choose columns 2, 4, 9 from the dataset (bonds + mkt + 5 industries)

data_ret = data[, c(2, 4, 9)] / 100
data_rf = data[, 3] / 100  # risk-free rate
data_excess_ret = sweep(data_ret, 1, data_rf, "-")   # Compute excess returns

# Compute mean vector and covariance matrix of excess returns (assuming they are constant)
Mu = colMeans(data_excess_ret)
Sigma = cov(data_excess_ret)

# Display the means and st.dev of the assets (in percentage) as well as the SR's
M = cbind(100*Mu, 100*sqrt(diag(Sigma)), Mu/sqrt(diag(Sigma)))
colnames(M) = c("Mean", "Stdev", "SR")
t(round(M, 3))

######################### 2. Portfolio analysis #########################
# Let us start with analyzing the properties of an equally weighted portfolio
nAssets = ncol(data_ret)
aNames = colnames(data_ret)
w = rep(1/nAssets, times = nAssets)
names(w) = aNames
t(round(w, 3))

# We will now compute the following of excess return: mean, std.dev, and Sharp
Stats = 100 * c(w %*% Mu, sqrt(w %*% Sigma %*% w ), (w %*% Mu) / sqrt(w %*% Sigma %*% w ) )
names(Stats) = c("Mean", "Stdev", "SR")
t(round(Stats, 3))


######################### 3. Tangency Portfolio #########################


##### Unconstrained #####
# We will compute the unconstrained tangency portfolio, i.e efficient portfolio of risky assets wih highest SR, allowed to take unlimited short position
# Most appropiate appraoch: Use mean vector and covariance matrix of excess return
w = uncTangent(Mu, Sigma)
t(round(w, ))

# Compute the mean, stdev and SR
uncT = list(w = w, mu=c(w %*% Mu), sd = c(sqrt(w %*% Sigma %*% w)),
            SR = c((w %*% Mu) / sqrt(w %*% Sigma %*% w)))

print(uncT)

##### Portfolio Constraints: No Short selling #####
w = conTangent(Mu, Sigma)
t(round(w,3 ))

# Compute the mean, stdev and SR
conT = list(w = w, mu = c(w %*% Mu), sd = c(sqrt(w %*% Sigma %*% w)), 
             SR = c((w %*% Mu)/(sqrt(w %*% Sigma %*% w))))
print(conT)

# Now we will specify portfolio constraints for each asset separately
lower = c(0.2, rep(0, times = nAssets - 1))
upper = c(1, rep(0.5, times = nAssets - 1))
w  = conTangent(Mu, Sigma, lower = lower, upper = upper)
t(round(w, 3))

# Compute the mean, stdev and SR
conT_alt = list(w = w, mu = c(w %*% Mu), sd = c(sqrt(w %*% Sigma %*% w)), 
            SR = c((w %*% Mu)/(sqrt(w %*% Sigma %*% w))))
print(conT_alt)



#############################################################################################################################################
############## Task 1 ##############

analyze_portfolio <- function(col_indices, data, rf_column) {
  # Extract the relevant columns and risk-free rate
  data_ret = data[, col_indices] / 100
  data_rf = data[, rf_column] / 100  # Risk-free rate
  
  # Compute excess returns
  data_excess_ret = sweep(data_ret, 1, data_rf, "-")
  
  # Compute mean vector and covariance matrix of excess returns
  Mu = colMeans(data_excess_ret)
  Sigma = cov(data_excess_ret)
  
  ##### Unconstrained #####
  w = uncTangent(Mu, Sigma)
  uncT = list(w = w, mu = c(w %*% Mu), sd = c(sqrt(w %*% Sigma %*% w)),
              SR = c((w %*% Mu) / sqrt(w %*% Sigma %*% w)))
  print("Unconstrained Portfolio:")
  print(uncT)
  
  ##### Constrained: No Short Selling #####
  w = conTangent(Mu, Sigma)
  conT = list(w = w, mu = c(w %*% Mu), sd = c(sqrt(w %*% Sigma %*% w)), 
              SR = c((w %*% Mu) / (sqrt(w %*% Sigma %*% w))))
  print("Constrained Portfolio (No Short Selling):")
  print(conT)
  
  # Portfolio constraints for each asset separately
  nAssets = length(col_indices)
  lower = c(0.2, rep(0, times = nAssets - 1))
  upper = c(1, rep(0.5, times = nAssets - 1))
  
  w = conTangent(Mu, Sigma, lower = lower, upper = upper)
  conT_alt = list(w = w, mu = c(w %*% Mu), sd = c(sqrt(w %*% Sigma %*% w)), 
                  SR = c((w %*% Mu) / (sqrt(w %*% Sigma %*% w))))
  print("Constrained Portfolio (Alternative Constraints):")
  print(conT_alt)
}

# Now we call the function with different column combinations
analyze_portfolio(c(2, 5, 8), data, 3)  
analyze_portfolio(c(4, 6, 7), data, 3)


#############################################################################################################################################
############## Task 2 ##############
# Find the optimal portfolio for an investor that does not want to short-sell and with risk-aversion 4.

# Parameters
gamma = 4
SR_tangent = conT$SR
sigma_tangent = conT$sd

# Leverage
l = SR_tangent / (gamma*sigma_tangent)
print(l)
portfolio_weights = l*uncT$w
final_weights = c(portfolio_weights, l-1)
print(final_weights)



######################### 4. No Risk-free rate #########################
# Below are functions that can be used to compute mean-variance optimal portfolio with or without constraints when there is no rf
# 1. for a given risk aversion  (gamma)
# 2. for a target mean
# 3. for a target st.dev as well as minimum variance portfolios

#minimum variance portfolios
w <- portMinVar(Sigma)                   #unconstrained
uncMV <- list(w = w, mu = c(w %*% Mu), sd = c(sqrt(w %*% Sigma %*% w)))
w <- portMinVar(Sigma, lower=0, upper=1) #no short-selling
conMV <- list(w = w, mu = c(w %*% Mu), sd = c(sqrt(w %*% Sigma %*% w)))

#no constraints: risk aversion, target mean, target st.dev.
w <- uncMeanVar_ra(Mu,Sigma,4)           #risk aversion 4
w <- uncMeanVar_mu(Mu,Sigma,0.01)        #target mean 1%
w <- uncMeanVar_sd(Mu,Sigma,0.05)         #target st.dev. 5%

#constraints:
w <- conMeanVar(Mu,Sigma, target="ra", t_val=4, lower=0, upper=1)
# w <- conMeanVar(Mu,Sigma,target="mu", t_val=0.008, lower=0, upper=1)
# w <- conMeanVar(Mu,Sigma,target="sd", t_val=0.04, lower=0, upper=1)


#############################################################################################################################################
############## Task 3 ##############

# Function to compute portfolio mean and stdev based on risk aversion
portfolio_stats <- function(gamma, Mu, Sigma) {
  w = uncMeanVar_ra(Mu, Sigma, gamma)  # Compute weights based on risk aversion
  mu_portfolio = c(w %*% Mu)           # Compute portfolio mean
  sigma_portfolio = c(sqrt(w %*% Sigma %*% w))  # Compute portfolio stdev
  return(list(mu = mu_portfolio, sd = sigma_portfolio))
}



# Target mean of 1% 
find_gamma_for_target_mean <- function(target_mu, Mu, Sigma) {
  target_function <- function(gamma) {
    stats = portfolio_stats(gamma, Mu, Sigma)
    return(stats$mu - target_mu)  # Solve for the difference between target mean and portfolio mean
  }
  gamma_solution = uniroot(target_function, interval = c(0.1, 10))$root  # Numerical method that iteratively find gamma 
                                          # which makes the difference between portfolio's mean and target equal to zero
  return(gamma_solution)
}

# Target standard deviation of 10%
find_gamma_for_target_stdev <- function(target_sd, Mu, Sigma) {
  target_function <- function(gamma) {
    stats = portfolio_stats(gamma, Mu, Sigma)
    return(stats$sd - target_sd)  # Solve for the difference between target stdev and portfolio stdev
  }
  gamma_solution = uniroot(target_function, interval = c(0.1, 10))$root  
  return(gamma_solution)
}

# Find the gamma for target mean of 1%
gamma_for_target_mean = find_gamma_for_target_mean(0.01, Mu, Sigma)
print(paste("Gamma for target mean of 1%:", round(gamma_for_target_mean, 4)))

# Find the gamma for target stdev of 10%
gamma_for_target_stdev = find_gamma_for_target_stdev(0.10, Mu, Sigma)
print(paste("Gamma for target stdev of 10%:", round(gamma_for_target_stdev, 4)))


######################### 5. Frontiers #########################
# Hypothetical rf'
rf = 0.002

# Range of expected excess returns
range_mu = seq(from=-max(Mu)*0.02, to=max(Mu)*1.2, by=max(Mu)/20)
n = length(range_mu)

# Function to compute the optimal portfolio for a given target mean "mu"
compute_portfolio <- function(mu, Mu, Sigma, rf, uncT) {
  # Without risk-free asset (wo rf)
  w = uncMeanVar_mu(Mu, Sigma, mu)
  
  # With risk-free asset (with rf)
  lev = mu / uncT$mu
  # Return portfolio statistics: Expected return, standard deviation, rf + mu, leverage-adjusted stdev
  return(c(rf + w %*% Mu, sqrt(w %*% Sigma %*% w), rf + mu, uncT$sd * lev))
}

# Use sapply to calculate results across the entire range of target excess returns
M = t(sapply(range_mu, compute_portfolio, Mu = Mu, Sigma = Sigma, rf = rf, uncT = uncT))

# plot assets
plot(sqrt(diag(Sigma)), rf+Mu, col="blue", xlab="St.Dev", ylab="Expc.ret",
     main="Efficient Frontiers", xlim=c(0,1.01*max(M[,2])), ylim=c(0, 1.01*max(M[,1])))

# Add efficient frontier without rf
lines(M[,2], M[, 1], col="black", lwd=2, lty=2)

# Add efficient frontier with rf
lines(M[,4], M[, 3], col="black", lwd=2, lty=1)

# Add tangency portfolio
points(uncT$sd, rf + uncT$mu, col="black", cex=2)

# Add minimum variance portfolio 
points(uncMV$sd, rf + uncMV$mu, col="purple", cex=2)

# Add grid
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

# MV efficient portfolio with constraints
x = conTangent(Mu - conMV$mu, Sigma, lower = 0, upper = 1)
mu_max = c(x %*% Mu)
range_mu = seq(from = conMV$mu * 1.01, to = mu_max, by = (mu_max - conMV$mu)/20)
n = length(range_mu)

# Function to compute the optimal portfolio for a given target mean "mu"
compute_optimal_portfolio <- function(mu, Mu, Sigma, rf, conT) {
  # Without risk-free asset (wo rf)
  w = uncMeanVar_mu(Mu, Sigma, mu)
  
  # With risk-free asset (with rf)
  lev = mu / uncT$mu
  # Return portfolio statistics: Expected return, standard deviation, rf + mu, leverage-adjusted stdev
  return(c(rf + w %*% Mu, sqrt(w %*% Sigma %*% w), rf + mu, uncT$sd * lev))
}

# Use sapply to calculate results across the entire range of target excess returns
N <- t(sapply(range_mu, compute_optimal_portfolio, Mu = Mu, Sigma = Sigma, conT = conT, rf = rf))

# Add efficient frontier without rf
lines(N[,2], N[,1], col = "red", lwd = 1, lty = 2)

# Add efficient frontier with rf
lines(N[,4], N[,3], col = "red", lwd = 1, lty = 1)

# Add tangency portfolio
points(conT$sd, rf + conT$mu, col = "red", cex = 1.5)

# Add minimum variance portfolio
points(conMV$sd, rf + conMV$mu, col = "green", cex = 1.5)
