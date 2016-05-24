#
# create.mixtures.R
#
# Copyright (c) 2013-2013 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified Mar, 2013
# first written Mar, 2013
# Contains: createMix, plot4Nmixture, createMix.test

# Create a mixture in X individuals per N 'M-modal' normals
# Default: 40 individuals per 4 conditions using a 2-model normal
# E is a vector of length N specifying what is going on in a condition (off, qtl, mean)
createMix <- function(X = 40, M = 2, E = c("qtl", "qtl", "mean", "off"), sd = 0.2, randomize = FALSE){
  N <- length(E)
  pheno <- NULL; geno <- NULL; cond <- NULL
  for(n in 1:N){ for(m in 1:M){
      mean <- m; sd <- sd
      if(E[n] == "off") { mean <- 1; sd <- sd/4 }
      if(E[n] == "mean"){ mean <- sum(1:M)/M }
      pheno <- c(pheno, rnorm(X/M, mean, sd))
      geno <- c(geno, rep(m,X/M)) # REAL genotype
      cond <- c(cond, rep(n,X/M)) # REAL condition
    }
  }
  m <- cbind(pheno, cond, geno)
  colnames(m) <- c("p","c","g")
  if(randomize) m <- m[sample(nrow(m)), ]
  invisible(m)
}

# Plot 4 conditions created by createMix as histograms (randomize should be set to false
plot4Nmixture <- function(mixture, X, N = 4, breaks=sqrt(length(mixture))){
  op <- par(mfrow=c(2,2))
  for(n in 0:(N-1)){ 
    s <- ((X*n)+1); e <- X*(n+1)
    hist(mixture[s:e], main=paste0("[",s,":",e,"]"), breaks=breaks)
  }
  op <- par(mfrow=c(1,1))
}

# Tests Create mix and the plot4Nmixture function
createMix.test <- function(){
  X = 100; M = 2 ; E = c("qtl", "off", "off", "off")

  m <- createMix(X, M, E, 0.15)
  plot4Nmixture(m[,"p"], X)
  hist(m[,"p"], breaks=sqrt(length(m[,"p"])))

  model <- lm(m[,"p"] ~ as.factor(m[,"c"]) + m[,"g"] + as.factor(m[,"c"]):m[,"g"])
  anova(model)
}

