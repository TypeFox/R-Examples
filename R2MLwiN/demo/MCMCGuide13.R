############################################################################
#     MLwiN MCMC Manual
#
# 13  Ordered Categorical Responses . . . . . . . . . . . . . . . . . . .181
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

# 13.1 A level chemistry dataset . . . . . . . . . . . . . . . . . . . . 181

library(R2MLwiN)
# MLwiN folder
mlwin <- getOption("MLwiN_path")
while (!file.access(mlwin, mode = 1) == 0) {
  cat("Please specify the root MLwiN folder or the full path to the MLwiN executable:\n")
  mlwin <- scan(what = character(0), sep = "\n")
  mlwin <- gsub("\\", "/", mlwin, fixed = TRUE)
}
options(MLwiN_path = mlwin)

# User's input if necessary

## Read alevchem data
data(alevchem, package = "R2MLwiN")

alevchem$gcseav <- double2singlePrecision(alevchem$gcse_tot/alevchem$gcse_no - 6)

hist(alevchem$gcseav, breaks = 20)

# 13.2 Normal response models . . . . . . . . . . . . . . . . . . . . . .184

(mymodel <- runMLwiN(a_point ~ 1 + (1 | pupil), estoptions = list(EstM = 1), data = alevchem))

(mymodel <- runMLwiN(a_point ~ 1 + gcseav + I(gcseav^2) + I(gcseav^3) + gender + (1 | pupil), estoptions = list(EstM = 1,
  resi.store = TRUE), data = alevchem))

predCurves(mymodel, xname = "gcseav", group = "genderfemale")

# 13.3 Ordered multinomial modelling . . . . . . . . . . . . . . . . . . 186

##IGLS
(mymodel <- runMLwiN(logit(a_point, cons, 6) ~ 1, D = "Ordered Multinomial", data = alevchem))

##MCMC
(mymodel <- runMLwiN(logit(a_point, cons, 6) ~ 1, D = "Ordered Multinomial", estoptions = list(EstM = 1), data = alevchem))

# 13.4 Adding predictor variables . . . . . . . . . . . . . . . . . . . .191
##MCMC
(mymodel <- runMLwiN(logit(a_point, cons, 6) ~ 1 + gcseav[1:5] + I(gcseav^2)[1:5] + I(gcseav^3)[1:5] + gender[1:5],
  D = "Ordered Multinomial", estoptions = list(EstM = 1), data = alevchem))

# 13.5 Multilevel ordered response modelling . . . . . . . . . . . . . . 192

# Note: Establishment codes on their own do not uniquely identify schools.
# Schools are instead uniquely identified by LEA code, establishment ID
# combination. Thus, here we generated a unique school ID.

alevchem$school <- as.numeric(factor(paste0(alevchem$lea, alevchem$estab)))

##MCMC
(mymodel <- runMLwiN(logit(a_point, cons, 6) ~ 1 + gcseav[1:5] + I(gcseav^2)[1:5] + gender[1:5] + (1[1:5] | school),
  D = "Ordered Multinomial", estoptions = list(EstM = 1), data = alevchem))

##MCMC
(mymodel <- runMLwiN(logit(a_point, cons, 6) ~ 1 + gcseav[1:5] + I(gcseav^2)[1:5] + gender[1:5] + (1[1:5] + gcseav[1:5] | school),
 D = "Ordered Multinomial", estoptions = list(EstM = 1), data = alevchem))
sixway(mymodel@chains[, "RP2_var_Intercept_12345", drop = FALSE], acf.maxlag = 300, "sigma2v6")

##Increases iterations to 50,000
(mymodel <- runMLwiN(logit(a_point, cons, 6) ~ 1 + gcseav[1:5] + I(gcseav^2)[1:5] + gender[1:5] + (1[1:5] + gcseav[1:5] | school),
 D = "Ordered Multinomial", estoptions = list(EstM = 1, mcmcMeth = list(iterations = 50000)), data = alevchem))
sixway(mymodel@chains[, "RP2_var_Intercept_12345", drop = FALSE], acf.maxlag = 300, "sigma2v6")

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .128





############################################################################
