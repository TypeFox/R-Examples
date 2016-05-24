### Create parameter sets for Generalized Inverse Gaussian

### gigLargeParam
chis <- c(0.1,0.2,0.5,1,2,5,10,20,50,100)
psis <- chis
lambdas <- c(-2,-1,-0.5,0,0.1,0.2,0.5,1,2,5,10)
maxrows <- length(lambdas)*length(chis)*length(psis)
gigLargeParam <- matrix(nrow = maxrows, ncol = 3)
rownum <- 1
for (i in 1:length(lambdas)){
  for (j in 1:length(chis)){
    for (k in 1:length(psis)){
      gigLargeParam[rownum,] <- c(chis[j],psis[k],lambdas[i])
      rownum <- rownum + 1
    }
  }
}

### gigSmallParam
chis <- c(0.1,0.5,2,10,50)
psis <- chis
lambdas <- c(-0.5,0,0.5,1,5)
maxrows <- length(lambdas)*length(chis)*length(psis)
gigSmallParam <- matrix(nrow = maxrows, ncol = 3)
rownum <- 1
for (i in 1:length(lambdas)){
  for (j in 1:length(chis)){
    for (k in 1:length(psis)){
      gigSmallParam[rownum,] <- c(chis[j],psis[k],lambdas[i])
      rownum <- rownum + 1
    }
  }
}

### save(gigLargeParam, gigSmallParam, file = "gigParam.rda")
