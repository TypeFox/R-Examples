################################################################################
# 
# StartingParams.R
# Version 1.1
# 30/01/2015
#
# Updates:
#   30/01/2015: Thomas model added
#
# Dataframes containing the starting parameters for the optimisation procedure
# fitting model parameters to observed coarse-scale occupancy data:
#   Nachman   Nachman model
#   PL        Power Law model
#   Logis     Logistic model
#   Poisson   Poisson model
#   NB        Negative binomial model
#   GNB       Generalised negative binomial model
#   INB       Improved negative binomial model
#   FNB       Finite negative binomial model
#   Thomas    Thomas model
#
################################################################################

### Nachman model
ParamsNachman <- data.frame("C" = 0.01, "z" = 0.01)

### Power Law model
ParamsPL <- data.frame("C" = 0.01, "z" = 0.01)

### Logistic model
ParamsLogis <- data.frame("C" = 0.01, "z" = 0.01)

### Poisson model
ParamsPoisson <- data.frame("lambda" = 1e-8)

### Negative binomial model
ParamsNB <- data.frame("C" = 0.01, "k" = 0.01)

### Generalised negative binomial model
ParamsGNB <- data.frame("C" = 0.00001, "z" = 1, "k" = 0.01)

### Improved negative binomial model
ParamsINB <- data.frame("C" = 1, "r" = 0.01, "b" = 0.1)

### Finite negative binomial model
ParamsFNB <- data.frame("W" = 10, "k" = 10)

### Thomas model
ParamsThomas <- data.frame("rho" = 1e-8, "mu" = 10,"sigma" = 1)
