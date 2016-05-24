############################################################################
#     MLwiN MCMC Manual
#
# 18  Multivariate Normal Response Models and Missing Data . . . . . . . 263
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

# 18.1 GCSE science data with complete records only . . . . . . . . . . .264

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

## Read gcsecomp1 data
data(gcsecomp1, package="R2MLwiN")
#summary(gcsecomp1)
#cor(gcsecomp1[, c("written", "csework")])

# 18.2 Fitting single level multivariate models . . . . . . . . . . . . .265

## IGLS
(mymodel <- runMLwiN(c(written, csework) ~ 1 + (1 | student), D = "Multivariate Normal", estoptions = list(sort.ignore = TRUE), 
  data = gcsecomp1))
## MCMC
(mymodel <- runMLwiN(c(written, csework) ~ 1 + (1 | student), D = "Multivariate Normal", estoptions = list(EstM = 1, 
  sort.ignore = TRUE), data = gcsecomp1))

# 18.3 Adding predictor variables . . . . . . . . . . . . . . . . . . . .270

(mymodel <- runMLwiN(c(written, csework) ~ 1 + female + (1 | student), D = "Multivariate Normal", estoptions = list(EstM = 1, 
  sort.ignore = TRUE), data = gcsecomp1))

# 18.4 A multilevel multivariate model . . . . . . . . . . . . . . . . . 271

## Store residual chain at level 3: school
(mymodel <- runMLwiN(c(written, csework) ~ 1 + female + (1 | school) + (1 | student), D = "Multivariate Normal", estoptions = list(EstM = 1, 
  resi.store = TRUE, resi.store.levs = 2), data = gcsecomp1))

resi <- mymodel@resi.chains$resi_lev2
label <- 1:ncol(resi)

## highlight
hipos <- rep(0, 6)
hipos[1] <- which(levels(as.factor(gcsecomp1$school)) == 68137)
hipos[2] <- which(levels(as.factor(gcsecomp1$school)) == 68201)
hipos[3] <- which(levels(as.factor(gcsecomp1$school)) == 68711)
hipos[4] <- which(levels(as.factor(gcsecomp1$school)) == 60427)
hipos[5] <- which(levels(as.factor(gcsecomp1$school)) == 22710)
hipos[6] <- which(levels(as.factor(gcsecomp1$school)) == 67105)

par(mfrow = c(2, 1))
## Select u0
resi0 <- resi[, label[which(label%%2 == 1)]]
resi0mean <- apply(resi0, 2, mean)
resi0sd <- apply(resi0, 2, sd)
rankno0 <- order(resi0mean)
resi0.lo <- resi0mean - 1.4 * resi0sd
resi0.hi <- resi0mean + 1.4 * resi0sd
caterpillar(y = resi0mean[rankno0], x = 1:length(resi0mean), qtlow = resi0.lo[rankno0], qtup = resi0.hi[rankno0], 
  ylim = c(-24, 21), ylab = "cons.written", xlab = "rank")
abline(h = 0, lty = "dotted")
for (i in 1:6) points(x = which(rankno0 == hipos[i]), y = resi0mean[rankno0[which(rankno0 == hipos[i])]], pch = 22, 
  bg = i + 1)

## Select u1
resi1 <- resi[, label[which(label%%2 == 0)]]
resi1mean <- apply(resi1, 2, mean)
resi1sd <- apply(resi1, 2, sd)
rankno1 <- order(resi1mean)
resi1.lo <- resi1mean - 1.4 * resi1sd
resi1.hi <- resi1mean + 1.4 * resi1sd
caterpillar(y = resi1mean[rankno1], x = 1:length(resi1mean), qtlow = resi1.lo[rankno1], qtup = resi1.hi[rankno1], 
  ylim = c(-24, 21), ylab = "cons.csework", xlab = "rank")
abline(h = 0, lty = "dotted")
for (i in 1:6) points(x = which(rankno1 == hipos[i]), y = resi1mean[rankno1[which(rankno1 == hipos[i])]], pch = 22, 
  bg = i + 1)

par(mfrow = c(1, 1))
plot(resi0mean, resi1mean, pch = 24, bg = "black", xlab = "cons.written", ylab = "cons.csework")
for (i in 1:6) points(x = resi0mean[rankno0[which(rankno0 == hipos[i])]], y = resi1mean[rankno1[which(rankno1 == hipos[i])]], 
  pch = 24, bg = i + 1)

# 18.5 GCSE science data with missing records . . . . . . . . . . . . . .275

## Read gcsemv1 data
data(gcsemv1, package = "R2MLwiN")

(mymodel <- runMLwiN(c(written, csework) ~ 1 + female + (1 | student), D = "Multivariate Normal", estoptions = list(EstM = 1, 
  sort.ignore = TRUE), data = gcsemv1))

(mymodel <- runMLwiN(c(written, csework) ~ 1 + female + (1 | school) + (1 | student), D = "Multivariate Normal", estoptions = list(EstM = 1, 
  mcmcMeth = list(dami = 2)), data = gcsemv1))

# 18.6 Imputation methods for missing data . . . . . . . . . . . . . . . 280

# 18.7 Hungarian science exam dataset . . . . . . . . . . . . . . . . . .281

## Read hungary1 data
data(hungary1, package = "R2MLwiN")
summary(hungary1)

(mymodel <- runMLwiN(c(es_core, biol_core, biol_r3, biol_r4, phys_core, phys_r2) ~ 1 + female + (1 | school) + (1 | student),
 D = "Multivariate Normal", data = hungary1))

(mymodel <- runMLwiN(c(es_core, biol_core, biol_r3, biol_r4, phys_core, phys_r2) ~ 1 + female + (1 | school) + (1 | student),
 D = "Multivariate Normal", estoptions = list(EstM = 1, mcmcMeth = list(dami = c(0, 1000, 2000, 3000, 4000, 
  5000))), data = hungary1))


if (!require(mitools)) install.packages("mitools")
library(mitools)
mi <- imputationList(mymodel@imputations)
head(mymodel@data)
with(mi, fun=head)

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .128





############################################################################
