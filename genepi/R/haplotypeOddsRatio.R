# this function calculates the haplotype disease risk adjusted for covariates
# formula: the logistic regression formula
# data:  name of the data frame containing the variables for analysis
# gvar:  names of genotype variables in data
# svar:  name of stratification variable (population stratification?)
# 
haplotypeOddsRatio <- function(formula, gtypevar, data, stratvar=NULL, nsim=100, tol=1e-8) {
  call <- match.call()
# remove missing
  data <- na.omit(data)
# number of subjects
  nsubjects <- nrow(data)
# number of loci
  k <- length(gtypevar)
# get the response variable
  disease <- model.response(model.frame(formula, data))
# setup the data and formula for analysis
  newdata <- data
  newformula <- update(formula, . ~ haplocopynum + .)
# genotype matrix indicating 0/1/2 copies of wt for each locus
  if(any(is.na(match(gtypevar, names(data))))) stop("one or more genotype variables not in the data frame")
  gtype <- as.matrix(2-data[gtypevar])
# create strata variable if missing
  if (missing(stratvar)) {
    strat <- rep(1,nsubjects)
  } else {
    if(is.na(match(stratvar, names(data)))) stop("stratification variable not in the data frame")
    strat <- data[[stratvar]]
  }
# Number of possible haplotypepairs
  two2k <- 2^k
  zzz <- genhaplopairs(k)
  hetmarkercount <- rowSums(1*(gtype==1))
# the pattern of heterozygous loci is converted from binary to a number
  g1idx <- (1*(gtype==1))%*%(2^(0:(k-1)))
# the pattern of wildtype homozygous loci is converted from binary to a number
  g2idx <- (1*(gtype==2))%*%(2^(0:(k-1)))
# compute the number of disease haplotypes the subject posses
  haplocopynum <- rep(0,nsubjects)
# if no locus is heterozygous and no locus is wt homozygous then 2 copies
  haplocopynum[hetmarkercount==0 & g2idx==0] <- 2
# if one locus is heterozygous and no locus is wt homozygous then 1 copy
  haplocopynum[hetmarkercount==1 & g2idx==0] <- 1
# if more than one heterozygous locus and no wt homozygous can be 0 or 1 copy
  newdata$haplocopynum <- as.factor(haplocopynum)
  ustrat <- unique(strat)
  copynumambig <- (hetmarkercount>1 & g2idx==0)
  namb <- sum(copynumambig)
  sortindex <- order(1*copynumambig, strat)
  copynumambig <- copynumambig[sortindex]
  strat <- strat[sortindex]
  disease <- disease[sortindex]
  g1idx <- g1idx[sortindex]
  g2idx <- g2idx[sortindex]
  newdata <- newdata[sortindex,]
  newdata1 <- newdata[copynumambig,]
  newdata1$haplocopynum <- rep(1,namb)
  newdata <- rbind(newdata, newdata1)
  p1copy <- rep(0,nsubjects)
  for(i in ustrat) {
#  for each population stratum get the controls
    i0 <- which(strat==i & disease==0)
    n0 <- length(i0)
#  compute the haplotype frequency under HWE  
    hf0 <- hwehaplofreq(n0, two2k, g1idx[i0], g2idx[i0], zzz$g1tbl, zzz$hpair)
#  for all subjects with ambiguous haplotypes in the stratum (includes cases)
    i1 <- which(strat==i & copynumambig)
    n1 <- length(i1)
#  compute the probability of 1 copy (case probabilities will be updated later)
    p1copy[i1] <- probonecopy(n1, hf0, g1idx[i1], zzz$g1tbl, zzz$hpair)
  }
# stop R CMD check warning: no visible binding for global variable ‘hfwts’
  hfwts <- NA
  newdata$hfwts <- c(1-p1copy,p1copy[copynumambig])
# haplotype ambiguous subjects will be coded as both copy num = 0 & 1
# the probability of 1 copy needs to be updated only for cases
# cases who will be coded as 0
  iamb0 <- which(disease==1 & copynumambig)
# cases who will be coded as 1
  iamb1 <- iamb0 + namb
  xyz <- glm(newformula, family=quasibinomial(), data=newdata, weights=hfwts)
# coefficients for intercept and haplocopynum=1
  haplo.lnor <- (xyz$coefficients)[1:2]
  or1copy1 <- 1
  or1copy0 <- exp(haplo.lnor[2])
  q1copy0 <- p1copy[iamb0]
  while (abs(or1copy0 - or1copy1) > tol) {
#  need exp(theta_1 + X beta) to update probability of 1 copy
    haplo.lp <- (xyz$linear.predictors - haplo.lnor[1])[iamb1]
    q1copy <- exp(haplo.lp)*q1copy0/(exp(haplo.lp)*q1copy0 + (1-q1copy0))
    newdata$hfwts[iamb1] <- q1copy
    newdata$hfwts[iamb0] <- 1-q1copy
    xyz <- glm(newformula, family=quasibinomial(), data=newdata, weights=hfwts)
    haplo.lnor <- (xyz$coefficients)[1:2]
    or1copy1 <- or1copy0
    or1copy0 <- exp(haplo.lnor[2])
  }
  aic <- deviance <- rep(0, nsim)
  orfit <- list()
  orfit$call <- call
  orfit$coef <- xyz$coefficients
  simcoef <- matrix(0,nsim,length(xyz$coefficients))
  orfit$var <- matrix(0,length(xyz$coefficients),length(xyz$coefficients))
  newdata1 <- newdata[1:nsubjects,]
  p <- xyz$rank
  p1 <- 1L:p
  for (i in 1:nsim) {
    newdata1$haplocopynum[nsubjects+1-(namb:1)] <- 1*(runif(namb) > newdata1$hfwts[nsubjects+1-(namb:1)])
    xyz <- glm(newformula, family=binomial(), data=newdata1)
    simcoef[i,] <- xyz$coefficients
# instead of summary(xyz)$cov.unscaled use only the relevant piece
    Qr <- xyz$qr
    covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    orfit$var <- orfit$var + covmat.unscaled
# collect the aic and deviance from each run of simulated haplotypes
    aic[i] <- xyz$aic
    deviance[i] <- xyz$deviance
  }
  dimnames(orfit$var) <- list(names(orfit$coef), names(orfit$coef))
  orfit$deviance <- mean(deviance)
  orfit$aic <- mean(aic)
  orfit$null.deviance <- xyz$null.deviance
  orfit$df.null <- xyz$df.null
  orfit$df.residual <- xyz$df.residual
  orfit$var <- orfit$var/nsim + var(simcoef)
  orfit$coef <- cbind(orfit$coef, sqrt(diag(orfit$var)), orfit$coef/sqrt(diag(orfit$var)), 2*pnorm(-abs(orfit$coef/sqrt(diag(orfit$var)))))
  colnames(orfit$coef) <- c("Estimate", "Std. Error", "Z-stat", "p.value")
  class(orfit) <- "haploOR"
  orfit
}


print.haploOR <- function(x, ...) {
  cat("\nCall:\n")
  dput(x$call)
  cat("\n")
  cat("Coefficients:\n")
  print(x$coef)
  cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null);   ",
      x$df.residual, "Residual\n")
  cat("Null Deviance:     ", x$null.deviance, 
      "\tResidual Deviance:", x$deviance, "\n")
}
