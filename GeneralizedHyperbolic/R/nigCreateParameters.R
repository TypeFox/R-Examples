### Create parameter sets for NIG
###
###
mus <- c(-1,0,1,2)
deltas <- c(1,2,5,10)
xis <- rep(c(0.1,0.3,0.5,0.7,0.9), 1:5)
chis <- c(0,-0.25,0.25,-0.45,0,0.45,-0.65,-0.3,0.3,0.65,
          -0.85,-0.4,0,0.4,0.85)
xiChis <- cbind(xis, chis)
lambdas <- -0.5


### nigLargeParam
maxrows <- length(mus)*length(deltas)*NROW(xiChis)*length(lambdas)
nigLargeParam <- matrix(nrow=maxrows,ncol=5)
rownum <- 1
for (i in 1:length(mus)){
  for (j in 1:length(deltas)){
    for (k in 1:NROW(xiChis)){
      for (l in 1:length(lambdas)){
        xi <- xiChis[k,1]
        chi <- xiChis[k,2]
        alpha <- (1 - xi^2)/(deltas[j]*xi*sqrt(xi^2 - chi^2))
        beta <- alpha*chi/xi
        param <- c(mus[i],deltas[j],alpha,beta,lambdas[l])
        nigLargeParam[rownum,] <- param
        rownum <- rownum + 1
      }
    }
  }
}
nigLargeParam <- nigLargeParam[,-5]

### nigSmallParam
mus <- c(0,1)
deltas <- c(1,5)
xis <- rep(c(0.1,0.3,0.7), c(1,2,4))
chis <- c(0,-0.25,0.25,-0.65,-0.3,0.3,0.65)
xiChis <- cbind(xis, chis)
lambdas <- -0.5

maxrows <- length(mus)*length(deltas)*NROW(xiChis)*length(lambdas)
nigSmallParam <- matrix(nrow=maxrows,ncol=5)
rownum <- 1
for (i in 1:length(mus)){
  for (j in 1:length(deltas)){
    for (k in 1:NROW(xiChis)){
      for (l in 1:length(lambdas)){
        xi <- xiChis[k,1]
        chi <- xiChis[k,2]
        alpha <- (1 - xi^2)/(deltas[j]*xi*sqrt(xi^2 - chi^2))
        beta <- alpha*chi/xi
        param <- c(mus[i],deltas[j],alpha,beta,lambdas[l])
        nigSmallParam[rownum,] <- param
        rownum <- rownum + 1
      }
    }
  }
}
nigSmallParam <- nigSmallParam[,-5]

### nigLargeShape
mus <- 0
deltas <- 1
xis <- rep(c(0.1,0.3,0.5,0.7,0.9), 1:5)
chis <- c(0,-0.25,0.25,-0.45,0,0.45,-0.65,-0.3,0.3,0.65,
          -0.85,-0.4,0,0.4,0.85)
xiChis <- cbind(xis, chis)
lambdas <- -0.5

maxrows <- length(mus)*length(deltas)*NROW(xiChis)*length(lambdas)
nigLargeShape <- matrix(nrow=maxrows,ncol=5)
rownum <- 1
for (i in 1:length(mus)){
  for (j in 1:length(deltas)){
    for (k in 1:NROW(xiChis)){
      for (l in 1:length(lambdas)){
        xi <- xiChis[k,1]
        chi <- xiChis[k,2]
        alpha <- (1 - xi^2)/(deltas[j]*xi*sqrt(xi^2 - chi^2))
        beta <- alpha*chi/xi
        param <- c(mus[i],deltas[j],alpha,beta,lambdas[l])
        nigLargeShape[rownum,] <- param
        rownum <- rownum + 1
      }
    }
  }
}
nigLargeShape <- nigLargeShape[,-5]

### nigSmallShape
mus <- 0
deltas <- 1
xis <- rep(c(0.1,0.3,0.7), c(1,2,4))
chis <- c(0,-0.25,0.25,-0.65,-0.3,0.3,0.65)
xiChis <- cbind(xis, chis)
lambdas <- -0.5

maxrows <- length(mus)*length(deltas)*NROW(xiChis)*length(lambdas)
nigSmallShape <- matrix(nrow=maxrows,ncol=5)
rownum <- 1
for (i in 1:length(mus)){
  for (j in 1:length(deltas)){
    for (k in 1:NROW(xiChis)){
      for (l in 1:length(lambdas)){
        xi <- xiChis[k,1]
        chi <- xiChis[k,2]
        alpha <- (1 - xi^2)/(deltas[j]*xi*sqrt(xi^2 - chi^2))
        beta <- alpha*chi/xi
        param <- c(mus[i],deltas[j],alpha,beta,lambdas[l])
        nigSmallShape[rownum,] <- param
        rownum <- rownum + 1
      }
    }
  }
}
nigSmallShape <- nigSmallShape[,-5]

## save(nigLargeParam, nigSmallParam, nigLargeShape, nigSmallShape,
##      file = "nigParam.rda")
