############################################################################
#     MLwiN MCMC Manual
#
# 4   Other Features of Variance Components Models . . . . . . . . . . . .45
#
#     Browne, W.J. (2009) MCMC Estimation in MLwiN, v2.13. Centre for
#     Multilevel Modelling, University of Bristol.
############################################################################
#     R script to replicate all analyses using R2MLwiN
#
#     Zhang, Z., Charlton, C., Parker, R, Leckie, G., and Browne, W.J.
#     Centre for Multilevel Modelling, 2012
#     http://www.bristol.ac.uk/cmm/software/R2MLwiN/
############################################################################

library(R2MLwiN)
# MLwiN folder
mlwin <- getOption("MLwiN_path")
while (!file.access(mlwin, mode = 1) == 0) {
  cat("Please specify the root MLwiN folder or the full path to the MLwiN executable:\n")
  mlwin <- scan(what = character(0), sep = "\n")
  mlwin <- gsub("\\", "/", mlwin, fixed = TRUE)
}
options(MLwiN_path = mlwin)

## Read tutorial data
data(tutorial, package = "R2MLwiN")

## IGLS
(mymodel1 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), data = tutorial))

## Gibbs
(mymodel2 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1), data = tutorial))

# 4.1 Metropolis Hastings (MH) sampling for the variance components model 46

# 4.2 Metropolis-Hastings settings . . . . . . . . . . . . . . . . . . . .47

# 4.3 Running the variance components with Metropolis Hastings . . . . . .48

## MH Adaptive with defaults
(mymodel3 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1, mcmcMeth = list(fixM = 2, 
  residM = 2, Lev1VarM = 2)), data = tutorial))

sixway(mymodel3@chains[, "FP_standlrt", drop = FALSE], "beta_1")

## MH Scale Factor = 5.8
(mymodel4 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1, mcmcMeth = list(fixM = 2, 
  residM = 2, Lev1VarM = 2, adaption = 0)), data = tutorial))

aa <- cbind(mymodel1@FP, mymodel2@FP, mymodel4@FP, mymodel3@FP)
bb <- cbind(mymodel1@RP, mymodel2@RP, mymodel4@RP, mymodel3@RP)
ctable <- round(rbind(aa, bb), 3)
colnames(ctable) <- c("IGLS", "Gibbs", "MH", "MH Adaptive")
print(ctable)
rm(list = c("mymodel1", "mymodel2", "mymodel3", "mymodel4"))
# 4.4 MH cycles per Gibbs iteration . . . . . . . . . . . . . . . . . . . 49

# 4.5 Block updating MH sampling . . . . . . . . . . . . . . . . . . . . .49
## MH Adaptive with defaults
(mymodel5 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1, mcmcMeth = list(fixM = 3, 
  residM = 2, Lev1VarM = 2, rate = 40)), data = tutorial))

estimates <- mymodel5@chains
par(mfrow = c(3, 2))
plot(4951:nrow(estimates), estimates[4951:nrow(estimates), "deviance"], xlab = "iteration", ylab = expression(paste("Est. of deviance")), 
  type = "l")
plot(4951:nrow(estimates), estimates[4951:nrow(estimates), "FP_Intercept"], xlab = "iteration", ylab = expression(paste("Est. of ", 
  beta[0])), type = "l")
plot(4951:nrow(estimates), estimates[4951:nrow(estimates), "FP_standlrt"], xlab = "iteration", ylab = expression(paste("Est. of ", 
  beta[1])), type = "l")
plot(4951:nrow(estimates), estimates[4951:nrow(estimates), "RP2_var_Intercept"], xlab = "iteration", ylab = expression(paste("Est. of ", 
  sigma[u0]^2)), type = "l")
plot(4951:nrow(estimates), estimates[4951:nrow(estimates), "RP1_var_Intercept"], xlab = "iteration", ylab = expression(paste("Est. of ", 
  sigma[e0]^2)), type = "l")
rm(mymodel5)

# 4.6 Residuals in MCMC . . . . . . . . . . . . . . . . . . . . . . . . . 51

(mymodel6 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1, resi.store = TRUE, 
  resi.store.levs = 2, mcmcMeth = list(iterations = 5001)), data = tutorial))

resi.chain2 <- mymodel6@resi.chains$resi_lev2
sixway(resi.chain2[, 1, drop = FALSE], name = "u0_1")

# 4.7 Comparing two schools . . . . . . . . . . . . . . . . . . . . . . . 54

dif <- resi.chain2[, 1] - resi.chain2[, 2]
sixway(dif, name = "dif")
prop <- (dif > 0)
mean(prop)

# 4.8 Calculating ranks of schools . . . . . . . . . . . . . . . . . . . .55

u0rank <- apply(resi.chain2, 1, rank)
u0rankmn <- apply(u0rank, 1, mean)
u0ranklo <- apply(u0rank, 1, function(x) quantile(x, 0.025))
u0rankmd <- apply(u0rank, 1, median)
u0rankhi <- apply(u0rank, 1, function(x) quantile(x, 0.975))

plot(1:65, u0rankmd, ylim = c(0.5, 65.5), pch = 15, xlab = "School", ylab = "Rank")
points(1:65, u0ranklo, pch = 24, bg = "grey")
points(1:65, u0rankhi, pch = 25, bg = "grey")
for (i in 1:65) lines(rep(i, 2), c(u0ranklo[i], u0rankhi[i]))

## common caterpillar plot

rankno <- order(u0rankmn)
plot(1:65, u0rankmn[rankno], ylim = c(0.5, 65.5), pch = 15, xlab = "School", ylab = "Rank")
points(1:65, u0ranklo[rankno], pch = 24, bg = "grey")
points(1:65, u0rankhi[rankno], pch = 25, bg = "grey")
for (i in 1:65) lines(rep(i, 2), c(u0ranklo[rankno[i]], u0rankhi[rankno[i]]))

# 4.9 Estimating a function of parameters . . . . . . . . . . . . . . . . 58
estimates <- mymodel6@chains
isc <- estimates[, "RP2_var_Intercept"]/(estimates[, "RP2_var_Intercept"] + estimates[, "RP1_var_Intercept"])
summary(isc)
sixway(isc, "isc")
rm(mymodel6)

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . . 60





############################################################################
