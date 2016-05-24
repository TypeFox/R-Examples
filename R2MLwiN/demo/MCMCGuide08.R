############################################################################
#     MLwiN MCMC Manual
#
# 8   Running a Simulation Study in MLwiN . . . . . . . . . . . . . . . . 97
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

# 8.1 JSP dataset simulation study . . . . . . . . . . . . . . . . . . . .97

# 8.2 Setting up the structure of the dataset . . . . . . . . . . . . . . 98

library(R2MLwiN)
# MLwiN folder
mlwin <- getOption("MLwiN_path")
while (!file.access(mlwin, mode = 1) == 0) {
  cat("Please specify the root MLwiN folder or the full path to the MLwiN executable:\n")
  mlwin <- scan(what = character(0), sep = "\n")
  mlwin <- gsub("\\", "/", mlwin, fixed = TRUE)
}
options(MLwiN_path = mlwin)

if (!require(doParallel)) install.packages("doParallel")
library(doParallel)

# User's input if necessary

# 8.3 Generating simulated datasets based on true values . . . . . . . . 102

# 8.4 Fitting the model to the simulated datasets . . . . . . . . . . . .106

pupil <- 1:108
school <- c(rep(1, 18), rep(2, 18), rep(3, 18), rep(4, 18), rep(5, 18), rep(6, 18))
cons <- rep(1, 108)

ns <- 100
IGLS_array <- MCMC_array <- array(, c(3, 2, ns))
MCMC_median <- data.frame(RP2_var_Intercept = rep(0, ns), RP1_var_Intercept = rep(0, ns))
CounterMCMC <- rep(0, 3)
Actual <- c(30, 10, 40)

simu <- function(i) {
    u_short <- rnorm(6, 0, sqrt(Actual[2]))
    u <- rep(u_short, each = 18, len = 108)
    e <- rnorm(108, 0, sqrt(Actual[3]))
    resp <- Actual[1] * cons + u + e
    indata <- data.frame(cbind(pupil, school, cons, resp))
    simModelIGLS <- runMLwiN(resp ~ 1 + (1 | school) + (1 | pupil), data = indata)
    simModelMCMC <- runMLwiN(resp ~ 1 + (1 | school) + (1 | pupil), estoptions = list(EstM = 1),
    data = indata)

    quantile25 <- c(quantile(simModelMCMC@chains[, "FP_Intercept"], 0.025),
                   quantile(simModelMCMC@chains[, "RP2_var_Intercept"], 0.025),
                   quantile(simModelMCMC@chains[, "RP1_var_Intercept"], 0.025))
    quantile975 <- c(quantile(simModelMCMC@chains[, "FP_Intercept"], 0.975),
                   quantile(simModelMCMC@chains[, "RP2_var_Intercept"], 0.975),
                   quantile(simModelMCMC@chains[, "RP1_var_Intercept"], 0.975))

    list(MCMCarray=cbind(coef(simModelMCMC), diag(vcov(simModelMCMC))),
         MCMCmed=c(median(simModelMCMC@chains[, "RP2_var_Intercept"]),
                   median(simModelMCMC@chains[, "RP1_var_Intercept"])),
         IGLSarray=cbind(coef(simModelIGLS), diag(vcov(simModelIGLS))),
         quantiles=cbind(quantile25, quantile975))
}

RNGkind("L'Ecuyer-CMRG")
cl <- makeCluster(detectCores(logical = FALSE))
registerDoParallel(cl)
set.seed(1)
r <- foreach(i=1:ns, .packages="R2MLwiN") %dopar% {
  simu(i)
}

stopCluster(cl)
unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister()

for(i in 1:ns){
  if (Actual[1] > r[[i]]$quantiles[1,1] & Actual[1] < r[[i]]$quantiles[1,2]) {
    CounterMCMC[1] <- CounterMCMC[1] + 1
  }
  if (Actual[2] > r[[i]]$quantiles[2,1] & Actual[2] < r[[i]]$quantiles[2,2]) {
    CounterMCMC[2] <- CounterMCMC[2] + 1
  }
  if (Actual[3] > r[[i]]$quantiles[3,1] & Actual[3] < r[[i]]$quantiles[3,2]) {
    CounterMCMC[3] <- CounterMCMC[3] + 1
  }
  MCMC_array[, , i] <- r[[i]]$MCMCarray
  MCMC_median[i, ] <- r[[i]]$MCMCmed
  IGLS_array[, , i] <- r[[i]]$IGLSarray
}

aa <- sapply(1:ns, function(x) na.omit(stack(as.data.frame(IGLS_array[, , x])))$values)
counterIGLS <- rep(0, 3)
for (i in 1:ns) {
  if (Actual[1] > aa[1, i] - 1.96 * sqrt(aa[4, i]) && Actual[1] < aa[1, i] + 1.96 * sqrt(aa[4, i]))
  {
    counterIGLS[1] <- counterIGLS[1] + 1
  }
  if (Actual[2] > aa[2, i] - 1.96 * sqrt(aa[5, i]) && Actual[2] < aa[2, i] + 1.96 * sqrt(aa[5, i]))
  {
    counterIGLS[2] <- counterIGLS[2] + 1
  }
  if (Actual[3] > aa[3, i] - 1.96 * sqrt(aa[6, i]) && Actual[3] < aa[3, i] + 1.96 * sqrt(aa[6, i]))
  {
    counterIGLS[3] <- counterIGLS[3] + 1
  }
}
Percent_interval_coverage <- (counterIGLS/ns) * 100
Mean_across_simus <- round(c(mean(aa[1, ]), mean(aa[2, ]), mean(aa[3, ])), 2)
Percent_bias <- round(-100 * (1 - Mean_across_simus/Actual), 2)
IGLS_results <- cbind(Mean_across_simus, Actual, Percent_bias, Percent_interval_coverage)
rownames(IGLS_results) <- c("beta0", "sigma2_u", "sigma2_e")
Percent_interval_coverage <- (CounterMCMC/ns) * 100
bb <- sapply(1:ns, function(x) na.omit(stack(as.data.frame(MCMC_array[, , x])))$values)
Mean_across_simus <- round(c(mean(bb[1, ]), mean(bb[2, ]), mean(bb[3, ])), 2)
Percent_bias <- round(-100 * (1 - Mean_across_simus/Actual), 2)
MCMC_results <- cbind(Mean_across_simus, Actual, Percent_bias, Percent_interval_coverage)
rownames(MCMC_results) <- c("beta0", "sigma2_u", "sigma2_e")

# 8.5 Analysing the simulation results . . . . . . . . . . . . . . . . . 109

cat("Simulation results using IGLS\n")
IGLS_results
cat("Simulation results using MCMC\n")
MCMC_results

# Investigating median estimates with Gamma(epsilon, epsilon) priors

Mean_across_simus <- round(c(mean(MCMC_median$RP2_var_Intercept),
mean(MCMC_median$RP1_var_Intercept)), 2)
Actual <- tail(Actual, -1)
Percent_bias <- round(-100 * (1 - Mean_across_simus/Actual), 2)
Percent_interval_coverage <- tail(Percent_interval_coverage, -1)
MCMC_results2 <- cbind(Mean_across_simus, Actual, Percent_bias, Percent_interval_coverage)
rownames(MCMC_results2) <- c("sigma2_u", "sigma2_e")
cat("Simulation results based on median MCMC estimates\n")
MCMC_results2

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . . 96





############################################################################
