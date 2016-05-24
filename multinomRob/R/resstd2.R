#
#  multinomRob
#
#  Walter R. Mebane, Jr.
#  University of Michigan
#  http://www-personal.umich.edu/~wmebane
#  <wmebane@umich.edu>
#
#  Jasjeet Singh Sekhon 
#  UC Berkeley
#  http://sekhon.polisci.berkeley.edu
#  sekhon@berkeley.edu
#
#

# probfunc: matrix of estimated probabilities
mnl.probfunc <-  function(Y, Ypos, Xarray, tvec) {
  nobs <- dim(Y)[1]
  ncats <- dim(Y)[2]
  eta <- matrix(0,nobs,ncats)
  for (j in 1:ncats) {
    useobs <- Ypos[,j];
    if (dim(tvec)[1] == 1) {
      eta[useobs,j] <- exp(Xarray[useobs,,j] * tvec[,j]);
    }
    else {
      eta[useobs,j] <- exp(Xarray[useobs,,j] %*% tvec[,j]);
    }
  }
  return( c(1/(eta %*% rep(1,ncats))) * eta )
}#end of mnl.probfunc

## residual.generator:  raw and median-centered ortho-standardized residuals
residual.generator <- function (tvec, Y, Ypos, X, m)
  {

    y.prob <- mnl.probfunc(Y, Ypos, X, tvec);

    if (all(Ypos)) {
      Sres.raw <- res.std(Y, m, y.prob);
      Sres <- Sres.raw - median(Sres.raw);
    }
    else {
      nobs <- dim(Y)[1];
      ncats <- dim(Y)[2];
      Sres.raw <- matrix(0, nobs, ncats-1);
      hasall <- apply(Ypos, 1, sum) == ncats;
      nobsall <- sum(hasall);
      if (nobsall > 0) {
        Yuse <- matrix(Y[hasall,], nobsall, ncats);  # in case nobsall == 1
        puse <- matrix(y.prob[hasall,], nobsall, ncats);
        Sres.raw[hasall,] <- res.std(Yuse, c(Yuse %*% rep(1,ncats)), puse);
      }
      hasless <- (1:nobs)[!hasall];
      Sres.use <- matrix(TRUE, nobs, ncats-1);
      for (i in hasless) {
        usecats <- Ypos[i,];
        nlesscats <- sum(usecats);
        ocats <- 1:(nlesscats-1);
        Yuse <- matrix(Y[i,usecats], 1, nlesscats);
        puse <- matrix(y.prob[i,usecats], 1, nlesscats);
        Sres.raw[i,ocats] <- res.std(Yuse, c(Yuse %*% rep(1,nlesscats)), puse);
        Sres.use[i,nlesscats:(ncats-1)] <- FALSE;
      }
      Sres <- Sres.raw - median(Sres.raw[Sres.use]);
      Sres[!Sres.use] <- 0;
    }

    return(list(Sres.raw=Sres.raw,Sres=Sres));
  }

## res.std: compute orthogonalized and standardized residuals
##  based on Tanabe and Sagae (1992 JRSSB)
res.std <- function(y,m,p, print.level=0)
  {
    obs    <- dim(y)[1]
    ncats  <- dim(y)[2]
    summat <- matrix(0,ncats,ncats)
    for (i in 1:(ncats-1)) summat[(i+1):ncats,i] <- 1;
    q <- p %*% summat
    ## d:  nonzero values in the diagonal matrix of the Cholesky decomposition
    d <- matrix(0,obs,ncats-1)
    for (i in 1:(ncats-1)) {
      if (i==1) d[,i] <- p[,i]*q[,i];
      if (i>1)  d[,i] <- p[,i]*q[,i]/q[,i-1];
    }
    ## r: raw residuals
    r <- y-m*p;
    summat <- matrix(0,ncats,ncats)
    for (i in 1:ncats) summat[i,i:ncats] <- 1;
    rsum <- r %*% summat
    ## Or:  orthogonalized residuals
    Or <- matrix(0,obs,ncats-1)
    for (i in 1:(ncats-1)) {
      if (i==1) Or[,i] <- r[,i];
      if (i>1)  Or[,i] <- r[,i] + rsum[,i-1]*p[,i]/q[,i-1];
    }
    ## Sr: standardized residuals
    Sr <- Or/sqrt(m*d);

    if (print.level > 0) {
      cat("res.std: Or\n");
      print(Or);
      cat("res.std: m\n");
      print(m);
      cat("res.std: d\n");
      print(d); 
    }
    return(Sr);
  }

kth.smallest  <- function(SortVector, obs, k)
  {
    if (obs < k) {
      print("ERROR! obs < k in kth.smallest\n");
      return(-1);
    }
    return(.Call("kthSmallest",
                 as.double(SortVector), as.integer(obs), as.integer(k),
                 PACKAGE="multinomRob"));    
  } #end of kth.smallest

