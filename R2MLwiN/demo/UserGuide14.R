############################################################################
#     MLwiN User Manual
#
# 14  Multivariate Response Models . . . . . . . . . . . . . . . . . . . 211
#
#     Rasbash, J., Steele, F., Browne, W. J. and Goldstein, H. (2012).
#     A User's Guide to MLwiN, v2.26. Centre for Multilevel Modelling,
#     University of Bristol.
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


# 14.1 Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . .211

data(gcsemv1, package = "R2MLwiN")
summary(gcsemv1)

# 14.2 Specifying a multivariate model . . . . . . . . . . . . . . . . . 212

# 14.3 Setting up the basic model . . . . . . . . . . . . . . . . . . . .214

(mymodel1 <- runMLwiN(c(written, csework) ~ 1 + (1 | student), D = "Multivariate Normal", estoptions = list(sort.ignore = TRUE), 
  data = gcsemv1))

(mymodel2 <- runMLwiN(c(written, csework) ~ 1 + female + (1 | school) + (1 | student), D = "Multivariate Normal", 
  data = gcsemv1))

mymodel2@RP["RP2_cov_Intercept_written_Intercept_csework"]/sqrt(mymodel2@RP["RP2_var_Intercept_written"] * mymodel2@RP["RP2_var_Intercept_csework"])

mymodel2@RP["RP1_cov_Intercept_written_Intercept_csework"]/sqrt(mymodel2@RP["RP1_var_Intercept_written"] * mymodel2@RP["RP1_var_Intercept_csework"])

# 14.4 A more elaborate model . . . . . . . . . . . . . . . . . . . . . .219

(mymodel3 <- runMLwiN(c(written, csework) ~ 1 + female + (1 + female | school) + (1 | student), D = "Multivariate Normal", 
  data = gcsemv1))

(mymodel4 <- runMLwiN(c(written, csework) ~ 1 + female + (1 + female[1] | school) + (1 | student), D = "Multivariate Normal", 
  estoptions = list(resi.store = TRUE), data = gcsemv1))

mymodel4@RP["RP2_cov_Intercept_written_Intercept_csework"]/sqrt(mymodel4@RP["RP2_var_Intercept_written"] * mymodel4@RP["RP2_var_Intercept_csework"])

mymodel4@RP["RP2_cov_Intercept_written_femaleFemale_1"]/sqrt(mymodel4@RP["RP2_var_Intercept_written"] * mymodel4@RP["RP2_var_femaleFemale_1"])

mymodel4@RP["RP2_cov_Intercept_csework_femaleFemale_1"]/sqrt(mymodel4@RP["RP2_var_Intercept_csework"] * mymodel4@RP["RP2_var_femaleFemale_1"])

u0 <- mymodel4@residual$lev_2_resi_est_Intercept.written
u1 <- mymodel4@residual$lev_2_resi_est_Intercept.csework
u2 <- mymodel4@residual$lev_2_resi_est_femaleFemale.1

plot(u0, u0, asp = 1)
plot(u0, u1, asp = 1)
plot(u0, u2, asp = 1)
plot(u1, u1, asp = 1)
plot(u1, u2, asp = 1)
plot(u2, u2, asp = 1)

# 14.5 Multivariate models for discrete responses . . . . . . . . . . . .222

data(tutorial, package = "R2MLwiN")

tutorial$binexam <- as.integer(tutorial$normexam > 0)
tutorial$binlrt <- as.integer(tutorial$standlrt > 0)

(mymodel5 <- runMLwiN(c(logit(binexam), logit(binlrt)) ~ 1, D = c("Mixed", "Binomial", "Binomial"), estoptions = list(sort.ignore = TRUE), 
  data = tutorial))


# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .224

############################################################################
