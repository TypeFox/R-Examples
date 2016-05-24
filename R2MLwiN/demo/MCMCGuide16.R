############################################################################
#     MLwiN MCMC Manual
#
# 16  Multiple Membership Models . . . . . . . . . . . . . . . . . . . . 231
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

# 16.1 Notation and weightings . . . . . . . . . . . . . . . . . . . . . 232

# 16.2 Office workers salary dataset . . . . . . . . . . . . . . . . . . 232

library(R2MLwiN)
# MLwiN folder
mlwin <- getOption("MLwiN_path")
while (!file.access(mlwin, mode = 1) == 0) {
  cat("Please specify the root MLwiN folder or the full path to the MLwiN executable:\n")
  mlwin <- scan(what = character(0), sep = "\n")
  mlwin <- gsub("\\", "/", mlwin, fixed = TRUE)
}
options(MLwiN_path = mlwin)

## Read wage1 data
data(wage1, package = "R2MLwiN")

summary(wage1)
hist(wage1$earnings)
hist(wage1$logearn, breaks = 20)

(mymodel <- runMLwiN(logearn ~ 1 + age_40 + numjobs + (1 | id), estoptions = list(EstM = 1), data = wage1))

(mymodel <- runMLwiN(logearn ~ 1 + age_40 + numjobs + sex + parttime + (1 | id), estoptions = list(EstM = 1), data = wage1))

vars <- cbind(as.numeric(wage1$parttime) - 1, as.numeric(wage1$sex) - 1, wage1$numjobs)
colnames(vars) <- c("parttime", "sex", "numjobs")
round(cor(vars), 4)

# 16.4 Fitting multiple membership models to the dataset . . . . . . . . 237

tabulate(wage1$numjobs)

(mymodel <- runMLwiN(logearn ~ 1 + age_40 + sex + parttime + (1 | company) + (1 | id), estoptions = list(EstM = 1), 
  data = wage1))

## Multiple membership
(mymodel <- runMLwiN(logearn ~ 1 + age_40 + sex + parttime + (1 | company) + (1 | id), estoptions = list(EstM = 1, 
  mm = list(list(mmvar = list("company", "company2", "company3", "company4"), weights = list("weight1", "weight2", 
    "weight3", "weight4")), NA), resi.store = TRUE, resi.store.levs = 2), data = wage1))

# 16.5 Residuals in multiple membership models . . . . . . . . . . . . . 240

lencateg <- length(unique(wage1$company))
resi0 <- mymodel@resi.chains$resi_lev2
resi0mean <- apply(resi0, 2, mean)
resi0sd <- apply(resi0, 2, sd)

rankno <- order(resi0mean)
resi0.lo <- resi0mean - 1.4 * resi0sd
resi0.hi <- resi0mean + 1.4 * resi0sd
caterpillar(y = resi0mean[rankno], x = 1:length(resi0mean), qtlow = resi0.lo[rankno], qtup = resi0.hi[rankno], ylim = c(-1, 
  1.3))
abline(h = 0, lty = "dotted")

aa <- qqnorm(resi0mean, plot.it = FALSE)
plot(x = aa$x[rankno], y = resi0mean[rankno], pch = 24, bg = "black", xlab = "nscore", ylab = "cons")
abline(h = 0, lty = "dotted")

wage1$companyno54 <- (wage1$company == 54) + (wage1$company2 == 54) + (wage1$company3 == 54) + (wage1$company4 == 
  54)
wage1$companyno67 <- (wage1$company == 67) + (wage1$company2 == 67) + (wage1$company3 == 67) + (wage1$company4 == 
  67)

## Multiple membership
(mymodel <- runMLwiN(logearn ~ 1 + age_40 + sex + parttime + companyno54 + companyno67 + (1 | company) + (1 | id), 
  estoptions = list(EstM = 1, mm = list(list(mmvar = list("company", "company2", "company3", "company4"), weights = list("weight1", 
    "weight2", "weight3", "weight4")), NA)), data = wage1))

#  16.6 Alternative weights for multiple membership models . . . . . . . .243


## New weights
(mymodel <- runMLwiN(logearn ~ 1 + age_40 + sex + parttime + companyno54 + companyno67 + (1 | company) + (1 | id), 
  estoptions = list(EstM = 1, mm = list(list(mmvar = list("company", "company2", "company3", "company4"), weights = list("weight1", 
    "weight2", "weight3", "weight4")), NA)), data = wage1))

# 16.7 Multiple membership multiple classification (MMMC) models . . . . 244

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .128





############################################################################
