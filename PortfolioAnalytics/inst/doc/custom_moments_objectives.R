### R code from vignette source 'custom_moments_objectives.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: custom_moments_objectives.Rnw:56-58
###################################################
library(PortfolioAnalytics)
library(DEoptim)


###################################################
### code chunk number 2: custom_moments_objectives.Rnw:63-72
###################################################
data(edhec)

# Use the first 4 columns in edhec for a returns object
R <- edhec[, 1:4]
colnames(R) <- c("CA", "CTAG", "DS", "EM")
head(R, 5)

# Get a character vector of the fund names
funds <- colnames(R)


###################################################
### code chunk number 3: custom_moments_objectives.Rnw:78-88
###################################################
# Construct initial portfolio with basic constraints.
init.portf <- portfolio.spec(assets=funds)
init.portf <- add.constraint(portfolio=init.portf, type="full_investment")
init.portf <- add.constraint(portfolio=init.portf, type="long_only")

# Portfolio with standard deviation as an objective
SD.portf <- add.objective(portfolio=init.portf, type="risk", name="StdDev")

# Portfolio with expected shortfall as an objective
ES.portf <- add.objective(portfolio=init.portf, type="risk", name="ES")


###################################################
### code chunk number 4: custom_moments_objectives.Rnw:92-97
###################################################
sd.moments <- set.portfolio.moments(R, SD.portf)
names(sd.moments)

es.moments <- set.portfolio.moments(R, ES.portf)
names(es.moments)


###################################################
### code chunk number 5: custom_moments_objectives.Rnw:114-121
###################################################
sigma.robust <- function(R){
  require(MASS)
  out <- list()
  set.seed(1234)
  out$sigma <- cov.rob(R, method="mcd")$cov
  return(out)
}


###################################################
### code chunk number 6: custom_moments_objectives.Rnw:125-129
###################################################
opt.sd <- optimize.portfolio(R, SD.portf, 
                             optimize_method="ROI", 
                             momentFUN="sigma.robust")
opt.sd


###################################################
### code chunk number 7: custom_moments_objectives.Rnw:133-138
###################################################
weights <- extractWeights(opt.sd)
sigma <- sigma.robust(R)$sigma

sqrt(t(weights) %*% sigma %*% weights)
extractObjectiveMeasures(opt.sd)$StdDev


###################################################
### code chunk number 8: custom_moments_objectives.Rnw:145-150
###################################################
pasd <- function(R, weights, sigma, N=36){
  R <- tail(R, N)
  tmp.sd <- sqrt(as.numeric(t(weights) %*% sigma %*% weights))
  sqrt(12) * tmp.sd
}


###################################################
### code chunk number 9: custom_moments_objectives.Rnw:168-177
###################################################
# Construct initial portfolio with basic constraints.
pasd.portf <- portfolio.spec(assets=funds)
pasd.portf <- add.constraint(portfolio=pasd.portf, type="full_investment")
pasd.portf <- add.constraint(portfolio=pasd.portf, type="long_only")

# Portfolio with pasd as an objective
# Note how we can specify N as an argument
pasd.portf <- add.objective(portfolio=pasd.portf, type="risk", name="pasd", 
                            arguments=list(N=48))


###################################################
### code chunk number 10: custom_moments_objectives.Rnw:182-187
###################################################
opt.pasd <- optimize.portfolio(R, pasd.portf, 
                               optimize_method="DEoptim", 
                               search_size=5000, trace=TRUE, traceDE=0,
                               momentFUN="sigma.robust")
opt.pasd


###################################################
### code chunk number 11: custom_moments_objectives.Rnw:200-211
###################################################
CRRA <- function(R, weights, lambda, sigma, m3, m4){
  weights <- matrix(weights, ncol=1)
  M2.w <- t(weights) %*% sigma %*% weights
  M3.w <- t(weights) %*% m3 %*% (weights %x% weights)
  M4.w <- t(weights) %*% m4 %*% (weights %x% weights %x% weights)
  term1 <- (1 / 2) * lambda * M2.w
  term2 <- (1 / 6) * lambda * (lambda + 1) * M3.w
  term3 <- (1 / 24) * lambda * (lambda + 1) * (lambda + 2) * M4.w
  out <- -term1 + term2 - term3
  out
}


###################################################
### code chunk number 12: custom_moments_objectives.Rnw:215-222
###################################################
crra.moments <- function(R, ...){
  out <- list()
  out$sigma <- cov(R)
  out$m3 <- PerformanceAnalytics:::M3.MM(R)
  out$m4 <- PerformanceAnalytics:::M4.MM(R)
  out
}


###################################################
### code chunk number 13: custom_moments_objectives.Rnw:226-237
###################################################
# Construct initial portfolio with basic constraints.
crra.portf <- portfolio.spec(assets=funds)
crra.portf <- add.constraint(portfolio=crra.portf, type="weight_sum", 
                             min_sum=0.99, max_sum=1.01)
crra.portf <- add.constraint(portfolio=crra.portf, type="box",
                             min=0.05, max=0.4)

# Portfolio with crra as an objective
# Note how we can specify lambda as an argument
crra.portf <- add.objective(portfolio=crra.portf, type="return", name="CRRA", 
                            arguments=list(lambda=10))


###################################################
### code chunk number 14: custom_moments_objectives.Rnw:240-244
###################################################
opt.crra <- optimize.portfolio(R, crra.portf, optimize_method="DEoptim",
                                 search_size=5000, trace=TRUE, traceDE=0,
                                 momentFUN="crra.moments")
opt.crra


