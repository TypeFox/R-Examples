#
#  multinomRob
#
#  Walter R. Mebane, Jr.
#  Cornell University
#  http://macht.arts.cornell.edu/wrm1/
#  wrm1@macht.arts.cornell.edu
#
#  Jasjeet Singh Sekhon 
#  UC Berkeley
#  http://sekhon.polisci.berkeley.edu
#  sekhon@berkeley.edu
#
#  $Id: multinomTanh.R,v 1.15 2005/09/27 08:04:06 wrm1 Exp $
#
#
# mGNtanh, multinomTanh and robustified.leverage
#
#

## mGNtanh:  Gauss-Newton tanh estimator, for overdispersed grouped multinomial GLM
##   Y:  matrix of (overdispersed and contaminated) multinomial counts
##   Ypos:  matrix indicating which in Y are counts (TRUE) and which are not (FALSE).
##   Xarray:  array of regressors,
##      dim(Xarray) = c(n observations, n parameters, n categories)
##   xvec:  vector to indicate all the coefficient parameters in the model
##      (parms by ncats):
##      It has a 1 for an estimated parameter, an integer >1 for an estimated
##      parameter constrained equal to another estimated parameter (all
##      parameters constrained to be equal to one another have the same integer
##      value in xvec) and a 0 otherwize.
##      example:
##      > xvec
##           [,1] [,2] [,3] [,4] [,5]
##      [1,]    1    1    1    1    0
##      [2,]    1    1    1    1    0
##      [3,]    1    1    1    1    0
##      [4,]    1    1    1    1    0
##   tvec: parms by ncats matrix
##   jacstack:  array of regressors,
##      dim(jacstack) = c(n observations, n UNIQUE parameters, n categories)
mGNtanh <- function(bstart, sigma2, resstart,
                      Y, Ypos, Xarray, xvec, tvec,
                      jacstack,itmax=100,print.level=0) {

  ## mcholeskyL:  lower-tri Cholesky matrix for multinomial;
  ##   p is vector of probs,
  mcholeskyL <- function(p) {
    n <- length(p);
    q <- rep(0,n);
    L <- diag(rep(1,n));
    for (i in 1:n) {
      if (i>1) L[i,1:(i-1)] <- -p[i]/q[1:(i-1)];
      if (i<n) q[i] <- 1 - sum(p[1:i]);
    }
    return(L);
  }

  ## mcholeskyLinv:  lower-tri inverse Cholesky matrix for multinomial;
  ##   p is vector of probs,
  mcholeskyLinv <- function(p) {
    n <- length(p);
    q <- rep(0,n);
    Linv <- diag(rep(1,n));
    for (i in 1:n) {
      if (i>1) Linv[i,1:(i-1)] <- p[i]/q[i-1];
      if (i<n) q[i] <- 1 - sum(p[1:i]);
    }
    return(Linv);
  }

  ## mcholeskyD:  Cholesky diagonal;
  ##   p is vector of probs,
  mcholeskyD <- function(p) {
    n <- length(p);
    q <- d <- rep(0,n);
    for (i in 1:n) q[i] <- 1 - sum(p[1:i]);
    d[1] <- p[1]*q[1];
    if (n>2) d[2:n] <- p[2:n]*q[2:n]/q[1:(n-1)];
    return(d);
  }

  ## probfunc: matrix of estimated probabilities
  probfunc <-
    function(Y, Ypos, Xarray, tvec) {
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
    }
  ## scorefunc:  score matrix
  scorefunc <-
    function(Ypos, nobs, nparms, phat, N, presmat, wmat, jacstack) {
      scoremat <- matrix(0,nparms,nobs);
      for (i in 1:nobs) {
        usecats <- Ypos[i,];
        nlesscats <- sum(usecats);
        pc <- mcholeskyL(phat[i,usecats]) ;
        pci <- mcholeskyLinv(phat[i,usecats]) ;
        adj <- t(pci) %*% diag(c(wmat[i,1:(nlesscats-1)],1)) %*% t(pc)
        ##     adj <- t(pci) %*% diag(c(wmat[i,1:(nlesscats-1)],1))
        if (dim(jacstack)[2]==1) {
          scoremat[,i] <- 
            N[i] * presmat[i,usecats] %*% adj %*% jacstack[i,,usecats] ;  ## weighted
        }
        else {
          scoremat[,i] <- 
            N[i] * presmat[i,usecats] %*% adj %*% t(jacstack[i,,usecats]) ;  ## weighted
        }
        ## if (dim(jacstack)[2]==1) {
        ##   scoremat[,i] <- 
        ##     N[i] * presmat[i,usecats] %*% jacstack[i,,usecats] ;  ## unweighted
        ## }
        ## else {
        ##   scoremat[,i] <- 
        ##     N[i] * presmat[i,usecats] %*% t(jacstack[i,,usecats]) ;  ## unweighted
        ## }
      }
      return( scoremat )
    }
  ## hessianfunc:  hessian matrix with simple weights
  hessianfunc <-
    function(Ypos, nobs, ncats, nparms, phat, N, wmat, jacstack) {
      H <- matrix(0,nparms,nparms)
      for (i in 1:nobs) {
        usecats <- Ypos[i,];
        nlesscats <- sum(usecats);
        pc <- mcholeskyL(phat[i,usecats]) ;
        pD <- mcholeskyD(phat[i,usecats]) ;
        wpvmat <- pc %*% diag(pD * c(wmat[i,1:(nlesscats-1)],1)) %*% t(pc) ;  ## matrix-weighted varmat
        ##     wpvmat <- diag(pD * c(wmat[i,1:(nlesscats-1)],1)) ;  ## matrix-weighted varmat
        ##     wpvmat <- diag(phat[i,usecats])-outer(phat[i,usecats],phat[i,usecats]);  ## unweighted
        H0 <- N[i] * wpvmat;
        if (dim(jacstack)[2]==1) {
          H <- H + jacstack[i,,usecats] %*% H0 %*% jacstack[i,,usecats]
        }
        else {
          H <- H + jacstack[i,,usecats] %*% H0 %*% t(jacstack[i,,usecats])
        }
      }
      return( H )
    }
  ## hessianfunc2:  mean hessian matrix with weights squared
  hessianfunc2 <-
    function(Ypos, nobs, ncats, nparms, phat, N, wmat, jacstack) {
      H <- matrix(0,nparms,nparms)
      for (i in 1:nobs) {
        usecats <- Ypos[i,];
        nlesscats <- sum(usecats);
        pc <- mcholeskyL(phat[i,usecats]) ;
        pD <- mcholeskyD(phat[i,usecats]) ;
        wpvmat <- pc %*% diag(pD * c(wmat[i,1:(nlesscats-1)]^2,1)) %*% t(pc) ;  ## matrix-weighted varmat
        ##     wpvmat <- diag(pD * c(wmat[i,1:(nlesscats-1)]^2,1)) ;  ## matrix-weighted varmat
        ##     wpvmat <- diag(phat[i,usecats])-outer(phat[i,usecats],phat[i,usecats]);  ## unweighted
        H0 <- N[i] * wpvmat;
        if (dim(jacstack)[2]==1) {
          H <- H + jacstack[i,,usecats] %*% H0 %*% jacstack[i,,usecats]
        }
        else {
          H <- H + jacstack[i,,usecats] %*% H0 %*% t(jacstack[i,,usecats])
        }
      }
      return( H / nobs)
    }
  ## resfunc:  orthogonalized and standardized (for multinomial covariance) resids
  resfunc <-
    function(Y, Ypos, Xarray, tvec) {
      if (all(Ypos)) {
        r <- res.std(Y, c(Y %*% rep(1,dim(Y)[2])), probfunc(Y, Ypos, Xarray, tvec));
      }
      else {
        nobs <- dim(Y)[1];
        ncats <- dim(Y)[2];
        r <- matrix(0, nobs, ncats-1);
        phat <- probfunc(Y, Ypos, Xarray, tvec);
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
  psifunc <-
    function(arg) {
      ## Hampel, Rousseeuw and Ronchetti 1981.  constants are from Table 2, p. 645
      ##                                                                          effic.
      ## c <- 3.0;    k <- 5.0;    A <- 0.680593;    B <- 0.769313;    d <- 1.470089;  ## 87%
      c <- 4.0;    k <- 5.0;    A <- 0.857044;    B <- 0.911135;    d <- 1.803134;  ## 97%
      return(
             ifelse(abs(arg)<d, arg,
                    ifelse(abs(arg)<c,
                           sqrt(A*(k-1))*tanh(sqrt((k-1)*B^2/A)*(c-abs(arg))/2)*sign(arg), 0))
             )
    }
  ## compute weight matrix
  weights <-
    function(Y, Ypos, Xarray, tvec, sigma2,ncats) {
      res <- resfunc(Y,Ypos,Xarray,tvec);
      sres <- res/sqrt(sigma2);
      ipsi <- psifunc(sres);
      wNA <- w <- ifelse(sres==0,1,ipsi/sres);
      if (!all(Ypos)) {  # set weights for nonexistent categories to zero
        npos <- apply(Ypos,1,sum);
        nobs <- dim(Y)[1];
        ncats <- dim(Y)[2];
        nwcats <- dim(w)[2];
        for (i in 1:nobs) {
          if (npos[i] < ncats) {
            w[i,npos[i]:nwcats] <- 0;
            wNA[i,npos[i]:nwcats] <- NA;
          }
        }
      }
      list(w=w,ipsi=ipsi,wNA=wNA)
    }
  ## check convergence
  converged <-
    function(bnew,bold) {
      return( sqrt(sum((bnew-bold)^2)) < 1e-6*(sqrt(sum(bold^2)) + 1e-4) )
    }
  converged2 <-
    function(snew,sold) {
      return( sqrt(sum((snew-sold)^2)) <= 1e-8*sqrt(sum(sold^2)) )
    }

  ## begin data computations
  Y[!Ypos] <- 0;  # ensure noncounts are set to zero, for convenience
  ncats <- dim(Y)[2]
  sres <- resstart/sqrt(sigma2);
  wmat <- ifelse(sres==0,1,psifunc(sres)/sres);
  bvec <- bstart;
  ## begin definition of variables used in GNstep that do not change over iterations
  nobs <- dim(Y)[1]
  ncats <- dim(Y)[2]
#  catidx <- 1:(ncats-1);
  tvars.total <- dim(Xarray)[2]
  tvunique <- dim(jacstack)[2]
  mvec <- c(Y %*% rep(1,dim(Y)[2]));
  propmat <- Y / mvec;  ## transform observed counts to proportions
  ## end definition of variables used in GNstep that do not change over iterations

  LogLik <-
    function(Y,Ypos,wmatS,ipmatS,mvecS) {
      LL <- 0;  LLu <- 0;
      for (i in 1:nobs) {
        usecats <- Ypos[i,];
        nlesscats <- sum(usecats);
        pc <- mcholeskyL(ipmatS[i,usecats]) ;
        pD <- mcholeskyD(ipmatS[i,usecats]) ;
        dnum <- diag( pc %*% diag(pD * c(wmat[i,1:(nlesscats-1)],0)) %*% t(pc) );
        pvec <- ipmatS[i,usecats];
        dden <- diag( diag(pvec)-outer(pvec,pvec) );
        adj <- dnum/dden;
        LL <- LL - sum(adj * Y[i,usecats] * log(pvec));  ## negative loglikelihood
        LLu <- LLu - sum(Y[i,usecats] * log(pvec));  ## negative loglikelihood
      }
      ##   print(paste("unweighted:",LLu,"; weighted:",LL))
      return(list(LL=LL,LLu=LLu));
    }

  ## Newton algorithm given data and a weight matrix (wmat)
  GNstep <-
    function(bvec,wmatGN,sigma2GN,itmaxGN=100,print.level=0) {
      tvec <- mnl.xvec.mapping(forward=FALSE,xvec,tvec,bvec, ncats,tvars.total);
      itersGN <- 0;
      ## patterned after the Gauss-Newton algorithm in Gallant 1987, 28-29
      for (iGN in 1:itmaxGN) {
        itersGN <- itersGN + 1;
        bprev2 <- bvec;
        
        ipmat <- probfunc(Y, Ypos, Xarray, tvec);
        presmat <- propmat - ipmat;
        loglik <- LogLik(Y,Ypos,wmatGN,ipmat,mvec)$LL;
        if (iGN==1 & (print.level > 1) )
          print(paste("mGNtanh: -loglik initial:",loglik));
        score <-
          scorefunc(Ypos, nobs, tvunique, ipmat, mvec, presmat, wmatGN, jacstack);
        hess2 <-
          hessianfunc2(Ypos, nobs, ncats, tvunique, ipmat, mvec, wmatGN, jacstack);
        ehess2 <- try(eigen(hess2, symmetric=TRUE, only.values=TRUE));
        if (length(ehess2) == 1 || any(ehess2$val==0)) {  # error in eigen() or singular
          posdef <- FALSE;
        }
        else {
          posdef <- all(ehess2$values > 0);
        }
        gradient <- (score %*% rep(1,nobs)) / sqrt(sigma2GN) / nobs ;
        if (!posdef) {
          convflag <- FALSE;
          break;  ## quit if Hessian is not positive definite
        }
        bdiff <- c(solve(hess2, tol=.Machine$double.eps) %*% gradient) ;  ## one Newton step
        ## print(bdiff);
        for (lambda in c(10:6)/10) {
          blambda <- bvec + lambda * bdiff;
          tlambda <- mnl.xvec.mapping(forward=FALSE,xvec,tvec,blambda, ncats,tvars.total);
          plambda <- probfunc(Y, Ypos, Xarray, tlambda);
          logliklambda <- LogLik(Y,Ypos,wmatGN,plambda,mvec)$LL;
          if (!is.na(logliklambda) && logliklambda < loglik) break;
        }
        if (is.na(logliklambda) || logliklambda > loglik) for (lambda in 2^(-(1:45))) {
          blambda <- bvec + lambda * bdiff;
          tlambda <- mnl.xvec.mapping(forward=FALSE,xvec,tvec,blambda, ncats,tvars.total);
          plambda <- probfunc(Y, Ypos, Xarray, tlambda);
          logliklambda <- LogLik(Y,Ypos,wmatGN,plambda,mvec)$LL;
          if (!is.na(logliklambda) && logliklambda < loglik) break;
        }
        if (is.na(logliklambda) | is.na(loglik)) {
          convflag <- FALSE;
          break;
        }
        if (logliklambda < loglik) bvec <- blambda;
        tvec <- mnl.xvec.mapping(forward=FALSE,xvec,tvec,bvec, ncats,tvars.total);
        if (print.level > 1)  {
          cat("mGNtanh: ibvec: "); print(bvec)
        }
        convflag <- converged(bvec,bprev2) & converged2(logliklambda,loglik);
        ##     convflag <- convflag & all(abs(gradient) < 1e-9);
        if (convflag) break;
      }
      if (!posdef & print.level >= 0) {
        print("mGNtanh: Hessian is not positive definite");
      }
      if (print.level > 1 & posdef) {
        print(paste("mGNtanh: -loglik final: ",logliklambda));
      }
      LL2 <- ifelse(posdef, LogLik(Y,Ypos,wmatGN,plambda,mvec), NA);
      if (print.level > 1 & posdef) {
        print(paste("mGNtanh: -loglik:  weighted,",LL2$LL,";  unweighted,",LL2$LLu));
        print("mGNtanh: gradient:");  print(c(gradient));
        print("mGNtanh: bvec:");  print(bvec);
      }
      information <-
        hessianfunc(Ypos, nobs, ncats, tvunique, ipmat, mvec, wmatGN, jacstack);
      ehess2 <- try(eigen(information, symmetric=TRUE, only.values=TRUE));
      if (length(ehess2) == 1 || any(ehess2$val==0)) {  # error in eigen() or singular hessian
        formation <- NA;
      }
      else {
        formation <- solve(information, tol=.Machine$double.eps);
      }
      return(
             list(coefficients=bvec, tvec=tvec, formation=formation, score=score,
                  LLvals=LL2, convflag=convflag, iters=itersGN, posdef=posdef) );
    }
  error <- 0;
  iters <- 0;
  for (i in 1:itmax) {
    iters <- iters + 1;
    bprev <- bvec;
    wprev <- wmat;

    tvec <- mnl.xvec.mapping(forward=FALSE,xvec,tvec,bvec, ncats,tvars.total);
    wres <- resfunc(Y,Ypos,Xarray,tvec)*wmat;
    wobs  <- sum(wmat);
    tanhsigma2  <- sum(wres^2)/(wobs-tvunique);

    ## grouped multinomial:  estimate using Newton algorithm
    GNlist <- GNstep(bvec,wmat,tanhsigma2,print.level);
    error <- ifelse(GNlist$posdef,0,32);  ## error == 32 if hessian not posdefinite
    if (print.level > 1)
      print(paste("mGNtanh: number of Newton iterations", GNlist$iters));

    bvec <- GNlist$coeff;
    wlist <- weights(Y,Ypos,Xarray,GNlist$tvec,sigma2,ncats);
    wmat <- wlist$w;
#    wmatNA <- wlist$wNA;
    if (error>0 || converged(bvec,bprev) & converged(wmat,wprev)) break;
  }

  opg <- GNlist$score %*% t(GNlist$score) ;
  obsformation <- GNlist$formation ;
  badformation <- length(obsformation)==1 && is.na(obsformation);
  if (badformation) {
    rcovmat <- obsformation <- NA * opg;
  }
  else {
    rcovmat <- obsformation %*% opg %*% obsformation;
  }

  if (print.level > 1) {  
    if (badformation) {
      print(paste("mGNtanh: hessian determinant:",NA));
    }
    else {
      print(paste("mGNtanh: hessian determinant:",
        det(solve(obsformation, tol=.Machine$double.eps))));
    }
    if (any(is.na(opg))) {
      print(paste("mGNtanh: OPG determinant:", NA));
    }
    else {
      print(paste("mGNtanh: OPG determinant:", det(opg)));
    }
    print(paste("mGNtanh: tanh sigma^2:", tanhsigma2));
  }

#  if (tanhsigma2 > sigma2) error <- error + 1;  ## error: tanh sigma2 > LQD sigma2
  if (error<32 && sum(wmat) < nobs*(ncats-1)/2) error <- error + 2;  ## error: wgts are too small

  ## table of returned error values (indicated values add to give total error)
  ## 0    no errors
  ## 1   tanh sigma2 > LQD sigma2
  ## 2   sum of weights < nobs*(ncats-1)/2
  ## 32  Hessian not positive definite in the final Newton step

  return(
         list(coefficients=GNlist$tvec, coeffvec=GNlist$coeff, dispersion=sigma2,
              w=wlist$wNA, psi=wlist$ipsi,
              A=opg/tanhsigma2, B=obsformation, covmat=rcovmat,
              iters=iters, error=error,
              GNlist=GNlist, tanhsigma2=tanhsigma2,
              Y=Y, Ypos=Ypos, probmat=probfunc(Y, Ypos, Xarray, GNlist$tvec),
              jacstack=jacstack, Xarray=Xarray)
         )
}

multinomTanh <- function (Y, Ypos, X, jacstack, xvec, tvec, pop, s2,
                     xvar.labels, choice.labels,
                     print.level=0) {

    nobs  <- dim(Y)[1];
    ncats <- dim(Y)[2];
    nvars <- dim(X)[2];
    nvars.unique <- sum(xvec == 1) + length(unique(xvec[xvec>1]));


    beta.vector <- vector(mode="numeric", length=nvars.unique);
    beta.vector <- mnl.xvec.mapping(forward=TRUE,xvec, tvec, beta.vector,
                                   ncats,nvars);

    residuals <- residual.generator(tvec,Y,Ypos,X,pop);

    mtanh <-
      mGNtanh(beta.vector, s2, residuals$Sres,
               Y, Ypos, X, xvec, as.matrix(tvec),
               jacstack,print.level=print.level);

    tvec <- mnl.xvec.mapping(forward=FALSE,xvec,tvec,mtanh$coeffvec,
                             ncats,nvars);

    se.vec <- sqrt(diag(mtanh$covmat));
    se  <- xvec;
    se  <- as.data.frame(mnl.xvec.mapping(forward=FALSE,xvec,se,se.vec,
                                          ncats,nvars));

    if (any(is.na(mtanh$A))) {
      se.opg.vec  <- rep(NA, dim(mtanh$A)[1]);
    }
    else {
      se.opg.vec  <- sqrt(diag(ginv(mtanh$A)));
    }
    se.opg  <- xvec;
    se.opg  <- as.data.frame(mnl.xvec.mapping(forward=FALSE,xvec,se.opg,se.opg.vec,
                                              ncats,nvars));

    if (any(is.na(mtanh$B))) {
      se.hes.vec  <- rep(NA, dim(mtanh$B)[1]);
    }
    else {
      se.hes.vec <- sqrt(diag(mtanh$B));
    }
    se.hes  <- xvec;
    se.hes  <- as.data.frame(mnl.xvec.mapping(forward=FALSE,xvec,se.hes,se.hes.vec,
                                              ncats,nvars));        

    row.names(se) <- xvar.labels;
    names(se)     <- choice.labels;
    mtanh$se      <- se;

    row.names(se.opg) <- xvar.labels;
    names(se.opg)     <- choice.labels;
    mtanh$se.opg      <- se.opg;

    row.names(se.hes) <- xvar.labels;
    names(se.hes)     <- choice.labels;
    mtanh$se.hes      <- se.hes;

    if (any(is.na(mtanh$w)) || any(is.na(jacstack))) {
      w.Hdiag <- Hdiag <- NA;
    }
    else {
      Hdiag <- robustified.leverage(tvec, Y, Ypos, X, pop, ifelse(mtanh$w>0,1,0),jacstack);
      w.Hdiag <- as.data.frame(matrix(c(as.vector(1:nobs),
                                 signif(mtanh$w), signif(Hdiag)),ncol=(ncats-1)+(ncats-1)+1));
      names(w.Hdiag) <-
        c("name",paste("weights:",choice.labels[1:ncats-1],sep=""),
          paste("Hdiag:",choice.labels[1:ncats-1],sep=""));
      #cat("mtanh: weights, Hdiag (by choices)\n");
    }

    sigma2 <- mtanh$disp;
    sigma  <- sqrt(sigma2);    

    cr <- fn.region.results(tvec, Y, Ypos, X, pop, sigma, Hdiag);

    mtanh$coef <- tvec;
    if (any(is.na(w.Hdiag))) {
      weights <- NA;
    }
    else {
      weights <- w.Hdiag[,c(1:ncats)];
      j  <- ncats
      Hdiag   <- w.Hdiag[,c(1,(j+1):(j+j-1))];
    }

    return( list(mtanh= mtanh,
                 weights=weights,
                 Hdiag=Hdiag,
                 cr   = cr,
                 tvec =tvec,
                 residuals= residual.generator(tvec,Y,Ypos,X,pop) ) );
  } #multinomTanh


robustified.leverage <- function (tvec, Y, Ypos, Xarray, m, Win,jacstack) {
#tvec <- mout$mtanh$coef
#Xarray <- X;
#m <- TotalY;
#W <- ifelse(mtanh$w>0,1,0);
#Win  <- mtanh$w

  obs  <- dim(jacstack)[1];
  W  <- cbind(Win,rep(1,obs));
  
  nobs  <- dim(Y)[1];
#  nvars <- dim(Xarray)[2];
  ncats <- dim(Xarray)[3];

  y.prob <- mnl.probfunc(Y, Ypos, Xarray, tvec);
  allYpos <- all(Ypos);
  if (!allYpos) {  # put probs for existing categories in y.prob[i,1:npos]
    npos <- apply(Ypos,1,sum);
    for (i in 1:nobs) {
      if (npos[i] < ncats) {
        y.prob[i,1:npos[i]] <- y.prob[i,Ypos[i,]];
        y.prob[i,(npos[i]+1):ncats] <- 0;
      }
    }
  }
  
  Hdiag <- matrix(0,nobs,ncats-1);
  summat <- matrix(0,ncats,ncats)
  for (i in 1:(ncats-1)) summat[(i+1):ncats,i] <- 1;
  p <- y.prob;
  q <- p %*% summat;

  #tanabe sagae '92
  L.gen <- function(p,q,ncats) {
    L <- matrix(0,ncats,ncats)
    for (j1 in 1:ncats) {
      for (j2 in 1:ncats) {
        if (j1 > j2)  { L[j1,j2]  <- -p[j1]/q[j2]; }
        if (j1 == j2) { L[j1,j2] <- 1; }
        if (j1 < j2)  { L[j1,j2] <- 0; }
      }
    }
    return(L)
  } #end of L.gen
  
  invL <- function(p,q,ncats) {
    Linv <- matrix(0,ncats,ncats)
    for (j1 in 1:ncats) {
      for (j2 in 1:ncats) {
        if (j1 > j2)  { Linv[j1,j2]  <- p[j1]/q[j1-1]; }
        if (j1 == j2) { Linv[j1,j2] <- 1; }
        if (j1 < j2)  { Linv[j1,j2] <- 0; }
      }
    }
    return(Linv)
  } #end of invL

# d:  nonzero values in the diagonal matrix of the Cholesky decomposition
  # note:  d[j,i]==0 if fewer than i+1 alternatives exist for obs j
  d <- matrix(0,obs,ncats-1)
  for (i in 1:(ncats-1)) {
    if (i==1) d[,i] <- p[,i]*q[,i];
    if (i>1)  d[,i] <- p[,i]*q[,i]/q[,i-1];
  }  

  Center  <- matrix(0,dim(jacstack)[2],dim(jacstack)[2]);
  for (i in 1:obs)
    {
      V <- matrix(0,nrow=ncats,ncol=ncats);
      V <- diag(c( ifelse(d[i,]>0, 1/sqrt(m[i]*d[i,]), 0) ,0));
      #USE THIS??? V <- V*ncats^2;

      L  <-  L.gen(p[i,],q[i,],ncats);
      
      if (!allYpos) {  # put jacstack for existing categories in smp[,1:npos]
        smp <- jacstack[i,,];
        if (npos[i] < ncats) {
          smp[,1:npos[i]] <- smp[,Ypos[i,]];
          smp[,(npos[i]+1):ncats] <- 0;
        }
        C0  <- smp %*% t(L) %*% V %*% diag(W[i,]) %*% V %*% t(L) %*% t(smp);
      }
      else {
        if (dim(jacstack)[2]==1) {
          C0  <- jacstack[i,,] %*% t(L) %*% V %*% diag(W[i,]) %*%
                  V %*% t(L) %*% jacstack[i,,];
        }
        else {
          C0  <- jacstack[i,,] %*% t(L) %*% V %*% diag(W[i,]) %*%
                   V %*% t(L) %*% t(jacstack[i,,]);
        }
      }
      Center  <- Center + C0;
    } #end of i

  invCenter  <- ginv(Center);

  for (i in 1:obs)
    {
      V <- matrix(0,nrow=ncats,ncol=ncats);
      V <- diag(c( ifelse(d[i,]>0, 1/sqrt(m[i]*d[i,]), 0) ,0));
      L  <-  L.gen(p[i,],q[i,],ncats);

      if (!allYpos) {  # put jacstack for existing categories in smp[,1:npos]
        smp <- jacstack[i,,];
        if (npos[i] < ncats) {
          smp[,1:npos[i]] <- smp[,Ypos[i,]];
          smp[,(npos[i]+1):ncats] <- 0;
        }
        Hdiag[i,]  <- diag(V %*% t(L) %*% t(smp) %*%
                        invCenter %*% smp %*% t(L) %*% V)[1:(ncats-1)];
      }
      else {
        if (dim(jacstack)[2]==1) {
          Hdiag[i,]  <-
            diag(V %*% t(L) %*% jacstack[i,,] %*%
              invCenter %*% jacstack[i,,] %*% t(L) %*% V)[1:(ncats-1)];
        }
        else {
          Hdiag[i,]  <-
            diag(V %*% t(L) %*% t(jacstack[i,,]) %*%
              invCenter %*% jacstack[i,,] %*% t(L) %*% V)[1:(ncats-1)];
        }
      }
    } #end of i
  
# put the negative forecasting variance adjustment values in Hdiag[w==0]
  for (j in 1:(ncats-1)) {
    sindx  <- Win[,j]==0;
    Hdiag[sindx,j] <- -Hdiag[sindx,j];
  }
  
  return(Hdiag);
}    #end of robustified.leverage

#calculate region results (orthogonal)
fn.region.results <- function (tmp.vec, Y, Ypos, X, TotalY, sigma, Hdiag) {
  y.prob      <- mnl.probfunc(Y,Ypos,X,tmp.vec);

  if (all(Ypos)) {
    Sres.raw <- res.std(Y, TotalY, y.prob);
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
      Sres.raw[hasall,] <- res.std(Yuse, TotalY[hasall], puse);
    }
    hasless <- (1:nobs)[!hasall];
    for (i in hasless) {
      usecats <- Ypos[i,];
      nlesscats <- sum(usecats);
      ocats <- 1:(nlesscats-1);
      Yuse <- matrix(Y[i,usecats], 1, nlesscats);
      puse <- matrix(y.prob[i,usecats], 1, nlesscats);
      Sres.raw[i,ocats] <- res.std(Yuse, TotalY[i], puse);
    }
  }

  Hmax <- .9;

  standard <- Sres.raw / sigma;
  student  <- Sres.raw / (sigma * sqrt(1-ifelse(Hdiag<Hmax,Hdiag,Hmax)));

  return(list(pred=y.prob,student=student,standard=standard));
} #end robustified leverage
