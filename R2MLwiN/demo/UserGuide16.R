############################################################################
#     MLwiN User Manual
#
# 16  An Introduction to Simulation Methods of Estimation . . . . . . . .241
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


# 16.1 An illustration of parameter estimation with Normally distributed . .
# . . .data . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .242

data(height, package = "R2MLwiN")
summary(height)

hist(height$height)

1 - pnorm((200 - mean(height$height))/sd(height$height))

heightsim1 <- function() {
  heightsim <- 175.35 + 10.002 * qnorm(runif(100))
  c(pmean = mean(heightsim), pvar = var(heightsim))
}

set.seed(1)

# Note: To obtain estimates as close as possible to the manual, increase the
# number of reps to 10000.

simdata1 <- as.data.frame(t(replicate(1000, heightsim1())))
simdata1$iteration <- 1:nrow(simdata1)

plot(simdata1$iteration, simdata1$pmean, type = "l")

plot(density(simdata1$pmean))

quantile(simdata1$pmean, c(0.025, 0.975))

plot(simdata1$iteration, simdata1$pvar, type = "l")

plot(density(simdata1$pvar))

quantile(simdata1$pvar, c(0.025, 0.975))

heightsim2 <- function(variable) {
  samp <- sample(variable, replace = TRUE)
  c(npmean = mean(samp), npvar = var(samp))
}

simdata2 <- as.data.frame(t(replicate(1000, heightsim2(height$height))))
simdata2$iteration <- 1:nrow(simdata2)

plot(simdata2$iteration, simdata2$npmean, type = "l")

plot(density(simdata2$npmean))

quantile(simdata2$npmean, c(0.025, 0.975))

plot(simdata2$iteration, simdata2$npvar, type = "l")

plot(density(simdata2$npvar))

quantile(simdata2$npvar, c(0.025, 0.975))

# 16.2 Generating random numbers in MLwiN . . . . . . . . . . . . . . . .249

female <- as.integer(runif(100) <= 0.6)

height2 <- (1 - female) * (175 + 10 * qnorm(runif(100))) + female * (160 + 8 * qnorm(runif(100)))

hist(height2[female == 0])
hist(height2[female == 1])

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .253

############################################################################
