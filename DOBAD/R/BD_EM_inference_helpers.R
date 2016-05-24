# "_helpers" not exported in package namespace


############################# the helpers/workers:

### Helper for getting Information Matrix for Fully observed BD proc
# (ie expectation of negative 2nd deriv of log likelihood)
# This _doesn't_ check that the expectation is calculated with the same
# parameters that are passed in. Parameters are assumed to be true.
getBDinform.full.SC.manual <- function(ENplus, ENminus, L, m){
  result <- matrix(data=NA, nrow=2, ncol=2);
  result[1,1] = ENplus / L^2;
  result[2,2] = ENminus / m^2
  result[1,2] <- result[2,1] <- 0;
  result;
}

#Louis '82 paper gives method for computing information
# This computes I_(X|Y) given true parameters
getBDinform.lost.SC.manual <- function(ENplus, ENminus, EHoldtime,
                             ENplusSq, ENminusSq, EHoldtimeSq,
                             ENplusNminus, ENplusHoldtime, ENminusHoldtime,
                             L, m, beta.immig, T){
  result <- matrix(data=NA, nrow=2, ncol=2);
  result[1,1] <- ENplusSq / L^2 - 2*ENplusHoldtime / L + EHoldtimeSq -
    2*(ENplus/L - EHoldtime)*T*beta.immig + T^2 *beta.immig^2;
  result[2,2] <- ENminusSq / m^2 - 2*ENminusHoldtime / m + EHoldtimeSq;
  result[1,2] <- result[2,1] <- ENplusNminus /(L*m) - ENminusHoldtime/ m -
    ENminus*T*beta.immig/m - ENplusHoldtime/L + EHoldtimeSq + EHoldtime*T*beta.immig;
  result;  
}

#NOTE: This doesn't calculate the last term in (3.2) in Louis ('82)
#since it is 0 for any thetahat satisfying the likelihood equation.
## Careful about T: make sure it corresponds with ENplus, etc
getBDinform.PO.SC.manual <- function(ENplus, ENminus, EHoldtime,
                             ENplusSq, ENminusSq, EHoldtimeSq,
                             ENplusNminus, ENplusHoldtime, ENminusHoldtime,
                             L, m, beta.immig, T){

#  if (){
#    BDloglikelihood.PO(partialDat=partialDat, L=L,m=m,nu=beta.immig*L);
#  }
  
  tmp1 <- getBDinform.full.SC.manual(ENplus, ENminus, L, m)
  tmp2 <- getBDinform.lost.SC.manual(ENplus, ENminus, EHoldtime,
                                  ENplusSq, ENminusSq, EHoldtimeSq,
                                  ENplusNminus, ENplusHoldtime, ENminusHoldtime,
                                  L, m, beta.immig, T);
  ##print(paste("inform.full", tmp1));
  ##print(paste("inform.lost", tmp2));
  tmp1-tmp2;
}

#Estimates E(fnc(Nt+)), E(fnc(Nt-)), and E(fnc(Rt)), conditionally.
# fnc should be from R to R, accepting a vector of Reals as argument
#Method = "cond" or "marg". referring to simulating condiitonally or
# marginally (and then dividing only by sims that fit the condigtional requriemtns)
# (IF this were to be coded , you wouldnt want to just simulate marginaly, but rather
# simulate each piece separately)
#getBDsummaryExpecs <- function(N=100, T, L, m, modelParams, data, n.fft=1024, method="cond",
#                               fnc){
getBDsummaryExpecs <- function(sims, fnc=function(x){x}){
  fullSummaries <- sapply(sims, BDsummaryStats, simplify=TRUE);
  N <- length(sims);
  results <- apply(fullSummaries, 1, function(x){mean(fnc(x));}) #fnc needs take vector arg
  names(results) <- c("Nplus", "Nminus", "Holdtime");
  results;
}

getBDsummaryProdExpecs <- function(sims, getsd=FALSE){
  fullSummaries <- sapply(sims, BDsummaryStats, simplify=TRUE);
  N <- length(sims);
  prods <- apply(fullSummaries, 2,
                  function(trip){c(trip[1]*trip[2], trip[1]*trip[3], trip[2]*trip[3]);}
                  );
  results <- apply(prods, 1, mean);
  names(results) <- c("NplusNminus", "NplusHoldtime", "NminusHoldtime");
  if (getsd){
    results <- matrix(c(results,
                        apply(prods, 1, sd)/sqrt(length(sims))
                        ), nrow=2, byrow=TRUE)
    return(results);
  }
  else
    return(results);
}





BDPOloglikeGradSqr.CTMC_PO_many <- function(partialDat, L,m, beta, n.fft=1024){
  ##require(numDeriv) ## in depends
  myLike <- function(theta){
    BDloglikelihood.PO.CTMC_PO_many(partialDat=partialDat,
                                    L=theta[1],m=theta[2],nu=theta[1]*beta, n.fft=n.fft);
  }
  mygrad <- matrix(grad(myLike, x=c(L,m)))
  mygrad %*% t(mygrad)
}

BDPOloglikeGradSqr.CTMC_PO_1 <- function(partialDat, L,m, beta, n.fft=1024){
  ##require(numDeriv) ## in depends
  myLike <- function(theta){
    BDloglikelihood.PO.CTMC_PO_1(partialDat=partialDat,
                                 L=theta[1],m=theta[2],nu=theta[1]*beta, n.fft=n.fft);
  }
  mygrad <- matrix(grad(myLike, x=c(L,m)))
  mygrad %*% t(mygrad)
}

###PLAN   currently using this to debug, Might want to keep it around later tho
#### So after done code it as getBDinform.PO.2, which uses the
#### extra gradient term to be safe.
#### [i.e.  result$I is the information]
getBDinform.PO.0term.CTMC_PO_many <- function(partialData, Lhat,Mhat,beta.immig,
                                        delta=1e-3, n=1024, r=4, prec.tol=1e-12,
                                        prec.fail.stop=1){
  myGetBDInform <- function(ctmc1){
    ## getBDinform.PO(partialData=ctmc1, Lhat=Lhat,Mhat=Mhat,
    ##                          beta.immig=beta.immig, delta=delta, n=n,
    ##                          r=r, prec.tol=prec.tol,
    ##                          prec.fail.stop=prec.fail.stop)
    ## ## renamed getBDinform.PO to getBDinform.PO.SC
    getBDinform.PO.SC(partialData=ctmc1, Lhat=Lhat,Mhat=Mhat,
                      beta.immig=beta.immig, delta=delta, n=n,
                      r=r, prec.tol=prec.tol,
                      prec.fail.stop=prec.fail.stop)

  }
  my0term <- function(ctmc1){
    BDPOloglikeGradSqr.CTMC_PO_1(partialDat=ctmc1, L=Lhat,m=Mhat,beta=beta.immig,n.fft=n)
  }
  ####dont know "efficient" way to get an array out ...
  ##sum(sapply(partialData@BDMCsPO, myGetBDInform));
  ##sapply(partialData@BDMCsPO, myGetBDInform);
  L <- length(partialData@BDMCsPO)
  informationsArray <- array(dim=c(2,2,L));
  zeroTerms <- array(dim=c(2,2,L));
  informationsArrayWZero <- array(dim=c(2,2,L));
  for (i in 1:L){
    informationsArray[,,i] <-
      myGetBDInform(partialData@BDMCsPO[[i]]);
    zeroTerms[,,i] <- my0term(partialData@BDMCsPO[[i]]);
    informationsArrayWZero[,,i] <- informationsArray[,,i]+zeroTerms[,,i];
    ##print(myGetBDInform(partialData@BDMCsPO[[i]]))
  }

  ##normally sum, but for debugging dont.
  ##apply(informationsArray,c(1,2),sum)
  
  list(informationsArray, zeroTerms, informationsArrayWZero,
       I.apprx=apply(informationsArray,c(1,2),sum),
       zeroTerm=apply(zeroTerms,c(1,2),sum),
       I=apply(informationsArrayWZero,c(1,2),sum))
}

getBDinform.PO.0term.CTMC_1 <- function(partialData, Lhat,Mhat,beta.immig,
                                        delta=1e-3, n=1024, r=4, prec.tol=1e-12,
                                        prec.fail.stop=1){
  myGetBDInform <- function(ctmc1){
    ## getBDinform.PO(partialData=ctmc1, Lhat=Lhat,Mhat=Mhat,
    ##                          beta.immig=beta.immig, delta=delta, n=n,
    ##                          r=r, prec.tol=prec.tol,
    ##                          prec.fail.stop=prec.fail.stop)

    ## ##  haven't tested but changing from getBDinform.PO to getBDinform.PO.SC
    getBDinform.PO.SC(partialData=ctmc1, Lhat=Lhat,Mhat=Mhat,
                   beta.immig=beta.immig, delta=delta, n=n,
                   r=r, prec.tol=prec.tol,
                   prec.fail.stop=prec.fail.stop)

  }
  my0term <- function(ctmc1){
    BDPOloglikeGradSqr.CTMC_PO_1(partialDat=ctmc1, L=Lhat,m=Mhat,beta=beta.immig,n.fft=n)
  }
  ####dont know "efficient" way to get an array out ...
  ##sum(sapply(partialData@BDMCsPO, myGetBDInform));
  ##sapply(partialData@BDMCsPO, myGetBDInform);
  informationsArray <- array(dim=c(2,2));
  zeroTerms <- array(dim=c(2,2));
  informationsArrayWZero <- array(dim=c(2,2));
  informationsArray <-  myGetBDInform(partialData);
  zeroTerms <- my0term(partialData);
  informationsArrayWZero <- informationsArray+zeroTerms;
  list(informationsArray, zeroTerms, informationsArrayWZero)
}




######
######  observed information via direct numeric calculation of hessian
###### (using numderiv package)


## rename, getBDinform.PO.SC.numeric
getBDinform.numeric <- function(partialData, Lhat,Mhat,beta.immig,
                                delta=1e-3,
                                r=4,
                                n.fft=1024){
  ##require(numDeriv) ## in Depends field of namepsace
  myLike <- function(theta){ ##ctmcpo1 or ctmcpomany
    BDloglikelihood.PO(partialDat=partialData,
                       L=theta[1],m=theta[2],nu=theta[1]*beta.immig,
                       n.fft=n.fft);
  }
##  Ihat <- -matrix(hessian(myLike,x=c(Lhat,Mhat), method.args=list(d=delta, r=r)))
    Ihat <- -hessian(myLike,x=c(Lhat,Mhat), method.args=list(d=delta, r=r))
  return(Ihat)
}

## BDPOloglikeGradSqr.CTMC_PO_1 <- function(partialDat, L,m, beta, n.fft=1024){
##   require(numDeriv)
##   myLike <- function(theta){
##     BDloglikelihood.PO.CTMC_PO_1(partialDat=partialDat,
##                                  L=theta[1],m=theta[2],nu=theta[1]*beta, n.fft=n.fft);
##   }
##   mygrad <- matrix(grad(myLike, x=c(L,m)))
##   mygrad %*% t(mygrad)
## }
