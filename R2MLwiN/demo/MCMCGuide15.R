############################################################################
#     MLwiN MCMC Manual
#
# 15  Cross Classified Models . . . . . . . . . . . . . . . . . . . . . .215
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

# 15.1 Classifications and levels . . . . . . . . . . . . . . . . . . . .216

# 15.2 Notation . . . . . . . . . . . . . . . . . . . . . . . . . . . . .217

# 15.3 The Fife educational dataset . . . . . . . . . . . . . . . . . . .217

library(R2MLwiN)
# MLwiN folder
mlwin <- getOption("MLwiN_path")
while (!file.access(mlwin, mode = 1) == 0) {
  cat("Please specify the root MLwiN folder or the full path to the MLwiN executable:\n")
  mlwin <- scan(what = character(0), sep = "\n")
  mlwin <- gsub("\\", "/", mlwin, fixed = TRUE)
}
options(MLwiN_path = mlwin)

## Read xc1 data
data(xc1, package = "R2MLwiN")

(mymodel <- runMLwiN(attain ~ 1 + (1 | sid) + (1 | pupil), estoptions = list(EstM = 1), data = xc1))

# 15.4 A Cross-classified model . . . . . . . . . . . . . . . . . . . . .220

(mymodel <- runMLwiN(attain ~ 1 + (1 | sid) + (1 | pid) + (1 | pupil), estoptions = list(xc = TRUE, EstM = 1, resi.store = TRUE, 
  resi.store.levs = c(2, 3)), data = xc1))

# 15.5 Residuals . . . . . . . . . . . . . . . . . . . . . . . . . . . . 223

lencateg <- length(unique(xc1$sid))
resi.chain0 <- mymodel@resi.chains$resi_lev3
residual0 <- apply(resi.chain0, 2, mean)
rankno <- order(residual0)
plot(x = 1:lencateg, y = residual0[rankno], pch = 24, bg = "black", xlab = "rank", ylab = "cons")
abline(h = 0, lty = "dotted")

## Common caterpillar
#lencateg <- length(unique(xc1[["sid"]]))
#resi.chain0 <- mymodel@resi.chains$resi_lev3
#u0rank <- apply(resi.chain0,1,rank)
#u0rankmn <- apply(u0rank, 1,mean)
#u0ranklo <- apply(u0rank, 1, function(x) quantile(x,.025))
#u0rankmd <- apply(u0rank, 1,median)
#u0rankhi <- apply(u0rank, 1, function(x) quantile(x,.975))
#rankno <- order(u0rankmn)
#caterpillar(y=u0rankmn[rankno],x=1:lencateg,qtlow=u0ranklo[rankno],qtup=u0rankhi[rankno]],ylim=c(0,20))

lencateg <- length(unique(xc1$pid))
resi.chain1 <- mymodel@resi.chains$resi_lev2
residual1 <- apply(resi.chain1, 2, mean)
rankno <- order(residual1)
plot(x = 1:length(residual1), y = residual1[rankno], pch = 24, bg = "black", xlab = "rank", ylab = "cons")
abline(h = 0, lty = "dotted")

## Common caterpillar
#lencateg <- length(unique(xc1[["pid"]]))
#resi.chain1 <- mymodel@resi.chains$resi_lev2
#u0rank <- apply(resi.chain1,1,rank)
#u0rankmn <- apply(u0rank, 1,mean)
#u0ranklo <- apply(u0rank, 1, function(x) quantile(x,.025))
#u0rankmd <- apply(u0rank, 1,median)
#u0rankhi <- apply(u0rank, 1, function(x) quantile(x,.975))
#rankno <- order(u0rankmn)
#caterpillar(y=u0rankmn[rankno],x=1:lencateg,qtlow=u0ranklo[rankno],qtup=u0rankhi[rankno],ylim=c(0,150))

# 15.6 Adding predictors to the model . . . . . . . . . . . . . . . . . .225

(mymodel <- runMLwiN(attain ~ 1 + vrq + (1 | sid) + (1 | pid) + (1 | pupil), estoptions = list(xc = TRUE, EstM = 1, 
  resi.store = TRUE, resi.store.levs = c(2, 3)), data = xc1))

(mymodel <- runMLwiN(attain ~ 1 + vrq + sc + fed + med + choice + (1 | sid) + (1 | pid) + (1 | pupil), estoptions = list(xc = TRUE, 
  EstM = 1, resi.store = TRUE, resi.store.levs = c(2, 3)), data = xc1))

lencateg <- length(unique(xc1$sid))
resi.chain0 <- mymodel@resi.chains$resi_lev3
residual0 <- apply(resi.chain0, 2, mean)
rankno <- order(residual0)
plot(x = 1:lencateg, y = residual0[rankno], pch = 24, bg = "black", xlab = "rank", ylab = "cons")
abline(h = 0, lty = "dotted")

xc1$school19 <- as.integer(xc1$sid == 19)

(mymodel <- runMLwiN(attain ~ 1 + vrq + sc + fed + med + choice + school19 + (0 | sid) + (1 | pid) + (1 | pupil), 
  estoptions = list(xc = TRUE, EstM = 1), data = xc1))

# 15.7 Current restrictions for cross-classified models . . . . . . . . .229

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .128





############################################################################
