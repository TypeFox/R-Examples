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
#  $Id: multinomMLE.R,v 1.11 2005/09/27 08:04:06 wrm1 Exp $
#
## multinomMLE:  maximum likelihood estimator for grouped multinomial GLM, with overdispersion
##   Y:  matrix of (overdispersed and contaminated) multinomial counts
##   Ypos:  matrix indicating which in Y are counts (TRUE) and which are not (FALSE).
##   Xarray:  array of regressors,
##      dim(Xarray) = c(n observations, n parameters, n categories)
##   xvec:  vector to indicate all the coefficient parameters in the model
##      (parms by ncats):
##      It has a 1 for an estimated parameter and a 0 otherwize.
##      example:
##      > xvec
##           [,1] [,2] [,3] [,4] [,5]
##      [1,]    1    1    1    1    0
##      [2,]    1    1    1    1    0
##      [3,]    1    1    1    1    0
##      [4,]    1    1    1    1    0
##   tvec: parms by ncats matrix (matrix with LQD estimates):
##      example:
##      > tvec
##                          Buchanan        Nader     Gore     Bush Other
##      int               -0.1641034    1.0735560 3.363641 4.151853     0
##      p(r,dg,d,r)96      2.4413780    0.3269827 3.207104 1.676676     0
##      p(cr,g,cd,cr)RV00 15.5333800 1149.4130000 2.039766 1.761392     0
##      pCuban            -6.4083750    0.3546630 1.966287 2.795598     0
##   jacstack:  array of regressors,
##      dim(jacstack) = c(n observations, n UNIQUE parameters, n categories)
multinomMLE <- function(Y, Ypos, Xarray, xvec, 
                        jacstack, itmax=100,
                        xvar.labels, choice.labels, MLEonly=FALSE,print.level=0) {
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
    function(Ypos, nobs, nparms, N, presmat, jacstack) {
      scoremat <- matrix(0,nparms,nobs);
      for (i in 1:nobs) {
        usecats <- Ypos[i,];
        if (dim(jacstack)[2]==1) {
          scoremat[,i] <- 
            N[i] * presmat[i,usecats] %*% jacstack[i,,usecats] ;  ## unweighted
        }
        else {
          scoremat[,i] <- 
            N[i] * presmat[i,usecats] %*% t(jacstack[i,,usecats]) ;  ## unweighted
        }
      }
      return( scoremat )
    }
  ## hessianfunc:  hessian matrix
  hessianfunc <-
    function(Ypos, nobs, nparms, phat, N, jacstack) {
      H <- matrix(0,nparms,nparms)
      for (i in 1:nobs) {
        usecats <- Ypos[i,];
        pvec <- phat[i,usecats];
        wpvmat <- diag(pvec)-outer(pvec,pvec);  ## unweighted
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
        for (i in hasless) {
          usecats <- Ypos[i,];
          nlesscats <- sum(usecats);
          ocats <- 1:(nlesscats-1);
          Yuse <- matrix(Y[i,usecats], 1, nlesscats);
          puse <- matrix(phat[i,usecats], 1, nlesscats);
          r[i,ocats] <- res.std(Yuse, c(Yuse %*% rep(1,nlesscats)), puse);
        }
      }
      return( r );
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
  tvec  <- xvec
  ## begin definition of variables used in GNstep that do not change over iterations
  nobs <- dim(Y)[1]
  ncats <- dim(Y)[2]
#  catidx <- 1:(ncats-1);
  tvars.total <- dim(Xarray)[2]
  tvunique <- dim(jacstack)[2]
  mvec <- c(Y %*% rep(1,dim(Y)[2]));
  propmat <- Y / mvec;  ## transform observed counts to proportions
  ## end definition of variables used in GNstep that do not change over iterations

  #starting values
  bvec <- rep(1,tvunique)

  if ((nobs-tvunique) <=0 & !MLEonly)
    {
      stop("too few observations to estimate overdispersed MNL.  You may want to set the 'MLEonly' option to true.")
    }

  LogLik <-
    function(Y,Ypos,ipmatS,mvecS) {
      LLu <- -sum(Y[Ypos] * log(ipmatS[Ypos]));  ## negative loglikelihood
      ##   print(paste("unweighted:",LLu))
      return(list(LLu=LLu));
    }

  ## Newton algorithm given data
  GNstep <-
    function(bvec,itmaxGN=100) {
      tvec <- mnl.xvec.mapping(forward=FALSE,xvec,tvec,bvec, ncats,tvars.total);
      itersGN <- 0;
      ## patterned after the Gauss-Newton algorithm in Gallant 1987, 28-29
      for (iGN in 1:itmaxGN) {
        itersGN <- itersGN + 1;
        bprev2 <- bvec;
        
        ipmat <- probfunc(Y, Ypos, Xarray, tvec);
        presmat <- propmat - ipmat;
        loglik <- LogLik(Y,Ypos,ipmat,mvec)$LLu;
        if (print.level > 32 & iGN==1) print(paste("multinomMLE: -loglik initial:",loglik));
        score <-
          scorefunc(Ypos, nobs, tvunique, mvec, presmat, jacstack);
        hess2 <- hessianfunc(Ypos, nobs, tvunique, ipmat, mvec, jacstack) / nobs;
        posdef <- all(eigen(hess2, symmetric=TRUE, only.values=TRUE)$values > 0);
        gradient <- (score %*% rep(1,nobs)) / nobs ;
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
          logliklambda <- LogLik(Y,Ypos,plambda,mvec)$LLu;
          if (!is.na(logliklambda) && logliklambda < loglik) break;
        }
        if (is.na(logliklambda) || logliklambda > loglik) for (lambda in 2^(-(1:45))) {
          blambda <- bvec + lambda * bdiff;
          tlambda <- mnl.xvec.mapping(forward=FALSE,xvec,tvec,blambda, ncats,tvars.total);
          plambda <- probfunc(Y, Ypos, Xarray, tlambda);
          logliklambda <- LogLik(Y,Ypos,plambda,mvec)$LLu;
          if (!is.na(logliklambda) && logliklambda < loglik) break;
        }
        if (logliklambda < loglik) bvec <- blambda;
        tvec <- mnl.xvec.mapping(forward=FALSE,xvec,tvec,bvec, ncats,tvars.total);
        if (print.level > 32)  {
          cat("multinomMLE: ibvec: "); print(bvec)
        }
        convflag <- converged(bvec,bprev2) & converged2(logliklambda,loglik);
        ##     convflag <- convflag & all(abs(gradient) < 1e-9);
        if (convflag) break;
      }
      if (!posdef & print.level >= 0) {
        print("multinomMLE: Hessian is not positive definite");
      }
      if (print.level > 32 & posdef) {
        print(paste("multinomMLE: -loglik final: ",logliklambda));
      }
      LL2 <- ifelse(posdef, LogLik(Y,Ypos,plambda,mvec), NA);
      if (print.level > 32 & posdef) {
        print(paste("multinomMLE: -loglik:  unweighted,",LL2$LLu));
        print("multinomMLE: gradient:");  print(c(gradient));
        print("multinomMLE: bvec:");  print(bvec);
      }
      information <- hessianfunc(Ypos, nobs, tvunique, ipmat, mvec, jacstack);
      if (all(eigen(information, symmetric=TRUE, only.values=TRUE)$values > 0)) {
        formation <- solve(information, tol=.Machine$double.eps);
      }
      else {
        formation <- NA;
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

    tvec <- mnl.xvec.mapping(forward=FALSE,xvec,tvec,bvec, ncats,tvars.total);
    if(MLEonly)
      {
        sigma2 <- 1
      } else {
        sigma2  <- sum(resfunc(Y,Ypos,Xarray,tvec)^2)/(nobs-tvunique);
      }

    ## grouped multinomial:  estimate using Newton algorithm
    GNlist <- GNstep(bvec);
    error <- ifelse(GNlist$posdef,0,32);  ## error == 32 if hessian not posdefinite
    if (print.level > 32)
      print(paste("multinomMLE: number of Newton iterations", GNlist$iters));

    bvec <- GNlist$coeff;
    if (converged(bvec,bprev)) break;
  }

  opg <- GNlist$score %*% t(GNlist$score) ;
  obsformation <- GNlist$formation ;
  rcovmat <- obsformation %*% opg %*% obsformation;

  if (print.level > 2) {  
    if (length(obsformation)==1 && obsformation==NA) {
      print(paste("multinomMLE: hessian determinant:",NA));
    }
    else {
      print(paste("multinomMLE: hessian determinant:",
        det(solve(obsformation, tol=.Machine$double.eps))));
    }
    print(paste("multinomMLE: OPG determinant:", det(opg)));
  }

  ## table of returned error values (indicated values add to give total error)
  ## 0    no errors
  ## 32  Hessian not positive definite in the final Newton step

  se.hes.vec <- sqrt(diag(obsformation));
  if(MLEonly)
    {
      se.rcov.vec <- rep(NA, length(diag(rcovmat)));
      se.opg.vec  <- rep(NA, length(diag(rcovmat)));      
    } else {
      se.rcov.vec     <- sqrt(diag(rcovmat));
      se.opg.vec  <- sqrt(diag(ginv(opg/sigma2)));
    }
  se.vec     <- se.hes.vec

  se.opg <- xvec;
  se.hes <- xvec;
  se.rcov <- xvec;
  se     <- xvec;    

  se.hes <- as.data.frame(mnl.xvec.mapping(forward=FALSE,xvec,se.hes,se.hes.vec,
                                           ncats,tvars.total));
  se <- as.data.frame(mnl.xvec.mapping(forward=FALSE,xvec,se,se.vec,
                                       ncats,tvars.total));
  if(MLEonly)
    {
      se.opg <- NA
      se.rcov <- NA
    } else {
      se.opg <- as.data.frame(mnl.xvec.mapping(forward=FALSE,xvec,se.opg,se.opg.vec,
                                               ncats,tvars.total));
      se.rcov <- as.data.frame(mnl.xvec.mapping(forward=FALSE,xvec,se.rcov,se.rcov,
                                                ncats,tvars.total));

      row.names(se.opg) <- xvar.labels;
      names(se.opg)     <- choice.labels;

      row.names(se.rcov) <- xvar.labels;
      names(se.rcov)     <- choice.labels;  
    }
  tvec  <- as.data.frame(tvec)

  row.names(tvec) <- xvar.labels;
  names(tvec)     <- choice.labels;   

  row.names(se) <- xvar.labels;
  names(se)     <- choice.labels;

  row.names(se.hes) <- xvar.labels;
  names(se.hes)     <- choice.labels;

  return(
         list(coefficients=GNlist$tvec, coeffvec=GNlist$coeff, dispersion=sigma2,
              se=se, se.opg=se.opg, se.hes=se.hes, se.sw=se.rcov,
              se.vec=se.vec, se.opg.vec=se.opg.vec,se.hes.vec=se.hes.vec,se.sw.vec=se.rcov.vec,
              A=opg/sigma2, B=obsformation, covmat=rcovmat,
              iters=iters, error=error,
              GNlist=GNlist, sigma2=sigma2,
              Y=Y, Ypos=Ypos, fitted.prob=probfunc(Y, Ypos, Xarray, GNlist$tvec),
              jacstack=jacstack)
         )
}

