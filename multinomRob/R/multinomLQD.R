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
#  $Id: multinomLQD.R,v 1.3 2005/09/23 03:12:42 wrm1 Exp $
#

##
## Multinomial lqd functions.
##
## MORE IN C CODE

## fit.multinomial.C.lqd2:  compute the LQD interquartile difference spread estimate
##
fit.multinomial.C.lqd2 <- function(foo,X,Y,Ypos,xvec,tvec,ncats,nvars,nvars.unique,obs,TotalY) {
  tmp.vec     <-   mnl.xvec.mapping(forward=FALSE,
                                    xvec,
                                    tvec,
                                    foo,
                                    ncats,
                                    nvars);  
  
  if (all(Ypos)) {
    return(.Call("original_fit_lqd2",
                 as.integer(obs),
                 as.integer(ncats),
                 as.integer(nvars),
                 as.integer(nvars.unique),
                 as.double(tmp.vec),
                 as.double(Y),
                 as.double(X),
                 as.double(TotalY),
                 PACKAGE="multinomRob")
           );
  }
  else {
    lqd2 <- function(r,nparms) {
      obs <- length(r)
      h <- ceiling((obs+nparms)/2);
      hidx <- (h*(h-1))/2;
      dif <- abs(outer(r,r,"-")[outer(1:obs,1:obs,">")])
    #  qrt <- sort(dif)[hidx]  * 2.21914446599;
      qrt  <-  kth.smallest(dif, length(dif), hidx)* 2.21914446599;  
      return(qrt)
    }

    Y[!Ypos] <- 0;
    sres <- resfunc.lqd2(Y, Ypos, X, tmp.vec);
    return( lqd2(c(sres[!is.na(sres)]),nvars.unique) );
  }
} #end of fit.multinomial.C.lqd2

#fit.multinomial.lqd2 <- function(foo,X,Y,xvec,tvec,ncats,nvars,nvars.unique) 
#{
#  tmp.vec     <- mnl.xvec.mapping(forward=FALSE,xvec,tmp.vec,foo,
#                                  ncats,nvars);
#  y.prob      <- mnl.probfunc(Y,X,tmp.vec);
#
#  Sres.raw    <- res.std(Y,TotalY,y.prob)
#
#  fit.lqd2 <- lqd2(Sres.raw,nvars.unique);
#  return(fit.lqd2);
#} #end of fit.multinomial.lqd2

## resfunc.lqd2:  orthogonalized and standardized (for multinomial covariance) resids
resfunc.lqd2 <- function(Y, Ypos, Xarray, tvec) {
  if (all(Ypos)) {
    r <- res.std(Y, c(Y %*% rep(1,dim(Y)[2])), mnl.probfunc(Y, Ypos, Xarray, tvec));
  }
  else {
    nobs <- dim(Y)[1];
    ncats <- dim(Y)[2];
    r <- matrix(NA, nobs, ncats-1);
    phat <- mnl.probfunc(Y, Ypos, Xarray, tvec);
    hasall <- apply(Ypos, 1, sum) == ncats;
    nobsall <- sum(hasall);
    if (nobsall > 0) {
      Yuse <- matrix(Y[hasall,], nobsall, ncats);  # in case nobsall == 1
      puse <- matrix(phat[hasall,], nobsall, ncats);
      r[hasall,] <- res.std(Yuse, c(Yuse %*% rep(1,ncats)), puse);
    }
    hasless <- (1:nobs)[!hasall];
    for (i in hasless) {  # orthostd resids go into r[i,1:(nlesscats-1)]
      usecats <- Ypos[i,];
      nlesscats <- sum(usecats);
      Yuse <- matrix(Y[i,usecats], 1, nlesscats);
      puse <- matrix(phat[i,usecats], 1, nlesscats);
      r[i,1:(nlesscats-1)] <- res.std(Yuse, c(Yuse %*% rep(1,nlesscats)), puse);
    }
  }
  return( r );
}
