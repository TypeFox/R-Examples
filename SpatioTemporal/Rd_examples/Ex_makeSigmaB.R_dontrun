##First create some random locations
x <- rnorm(5)
y <- rnorm(5)

##compute distance matrix
D <- crossDist( cbind(x,y) )

##create a block diagonal matrix exponential covariance matrix
##with different range, sill, and nugget
pars <- list(c(.3,2), c(2,1), c(1,3))
nugget <- c(.5,0,1)

Sigma1 <- makeSigmaB(pars, D, type="exp", nugget=nugget)
##or using different covariance functions for each block
Sigma2 <- makeSigmaB(pars, D, type=c("exp","exp2","cubic"),
                     nugget=nugget)

##make a cross-covariance matrix
Dcross <- D[1:3,c(1,1,2,2)]
Sigma.cross <- makeSigmaB(pars, Dcross, type="exp", nugget=nugget,
                          ind2.to.1=c(1,1,2,2))

\dontshow{
  Sigma.alt <- matrix(0, length(pars)*dim(D)[1], length(pars)*dim(D)[1])
  for(i in 1:length(pars)){
    Ind <- (1:dim(D)[1]) + (i-1)*dim(D)[1]
    Sigma.alt[Ind, Ind] <- pars[[i]][2]*exp(-D/pars[[i]][1])
    diag(Sigma.alt[Ind, Ind]) <- diag(Sigma.alt[Ind, Ind])+nugget[i]
  }
  if( abs(max(Sigma1-Sigma.alt)) > 1e-13){
    stop("makeSigmaB: Results not equal, covariance")
  }
  Ind <- c(1,1,2,2)
  Sigma.alt.cross <- Sigma.alt[c(1:3,6:8,11:13),c(Ind, 5+Ind, 10+Ind)]
  if( abs(max(Sigma.cross-Sigma.alt.cross)) > 1e-13){
    stop("makeSigmaB: Results not equal, cross-covariance")
  }

  Sigma2.s <- makeSigmaB(pars, D, type=c("exp","exp2","cubic"),
                         nugget=nugget, sparse=TRUE)
  if( max(abs(Sigma2.s-Sigma2)) > 1e-14 ){
    stop("makeSigmaB: Sparse matrix not equal, covariance")
  }
  Sigma.cross.s <- makeSigmaB(pars, Dcross, type="exp", nugget=nugget,
                              ind2.to.1=c(1,1,2,2), sparse=TRUE)
  if( max(abs(Sigma.cross.s-Sigma.cross)) > 1e-14 ){
    stop("makeSigmaB: Sparse matrix not equal, cross-covariance")
  }
}

