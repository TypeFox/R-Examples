############################################################################
#     MLwiN User Manual
#
# 3   Residuals . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 37
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


# 3.1 What are multilevel residuals? . . . . . . . . . . . . . . . . . . .37

data(tutorial, package = "R2MLwiN")

(mymodel1 <- runMLwiN(normexam ~ 1 + (1 | school) + (1 | student), data = tutorial, estoptions = list(resi.store = TRUE)))



# 3.2 Calculating residuals in MLwiN . . . . . . . . . . . . . . . . . . .40

residuals <- mymodel1@residual$lev_2_resi_est_Intercept
residualsCI <- 1.96 * sqrt(mymodel1@residual$lev_2_resi_var_Intercept)
residualsRank <- rank(residuals)
rankno <- order(residualsRank)

caterpillar(y = residuals[rankno], x = 1:65, qtlow = (residuals - residualsCI)[rankno],
           qtup = (residuals + residualsCI)[rankno], xlab = 'Rank', ylab = 'Intercept')

# 3.3 Normal plots . . . . . . . . . . . . . . . . . . . . . . . . . . . .43

e0 <- mymodel1@residual$lev_1_resi_est_Intercept

e0std <- (e0 - mean(e0))/sd(e0)

e0rank <- rank(e0)

e0uniform <- (e0rank - 0.5)/length(e0rank)

e0nscore <- qnorm(e0uniform)

plot(e0nscore, e0std, asp = 1)


u0 <- mymodel1@residual$lev_2_resi_est_Intercept

u0std <- (u0 - mean(u0))/sd(u0)

u0rank <- rank(u0)

u0uniform <- (u0rank - 0.5)/length(u0rank)

u0nscore <- qnorm(u0uniform)

plot(u0nscore, u0std, asp = 1)

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . . 45


############################################################################
