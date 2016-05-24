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
#  $Id: multinomRob.R,v 1.16 2005/09/27 08:04:06 wrm1 Exp $


######################################################
## multinomRob: genoud and Gauss-Newton estimation  ##
######################################################

multinomRob <-
  function(model,
           data,
           starting.values=NULL,
           equality=NULL,  # list of lists of parameter equality constraints
           genoud.parms   = NULL, ## can be NULL, single, or list
           print.level    = 0,
           iter = FALSE,   # should we iterate
           maxiter = 10,  # maximum number of iterations before we stop
           multinom.t=1,  # 0=no, 1=yes, 2=force should we do multinom-t for starting values
           multinom.t.df=NA, #if set, the multivariate-t function is FORCED to use this DF
           MLEonly=FALSE
           ) 
  {

    #check input
    if (print.level < 0)
      print.level  <- 0;

    if(MLEonly!=0 & MLEonly!=1)
      MLEonly <- FALSE

    if (iter!=TRUE & iter!=FALSE)
      {
        stop("multinomRob(): illegal input for variable iter")
      }    
    if (maxiter < 0)
      {
        stop("multinomRob(): illegal input for variable maxiter")
      }
    if (multinom.t!=0 & multinom.t!=1 & multinom.t!=2)
      {
        stop("multinomRob(): illegal input for variable multinom.t")
      }
    if (!is.na(multinom.t.df) & multinom.t.df <= 0)
      {
        stop("multinomRob(): illegal input for variable multinom.t.df")
      }    

    converged.test <- function(old,new,tol=1e-8) {
      abs((old-new)/old) < tol;
    }
    genoud.fun <- function(z) {
      fit.multinomial.C.lqd2(z,X,Y,Ypos,xvec,tvec,ncats,nvars,nvars.unique,obs,TotalY)
    }

    ## ----------------- BEGIN: PARSING / MODEL BUILDING AUTOMATION ------------------

    ## ---- DATA X, Y and Ypos ----
    data.all <- get.xy(model, data, print.level=print.level);
    Y <- data.all$Y;
    ncats  <- dim(Y)[2];
    obs    <- dim(Y)[1];
    Ypos <- data.all$ypos;  # element is TRUE if y>=0, FALSE if y<0
    Y[!Ypos] <- 0;  # set count values to be ignored to zero, for convenience later
    X <- data.all$X;
    nvars <- dim(X)[2];

    ## ---- LABELS AND CONTROLS/FLAGS ----
    xvec <- matrix(0,nrow=nvars,ncol=ncats)
    colnames(xvec) <-  choice.labels  <- data.all$ynames;

    xvar.labels <- rep( "", nvars)
    for (i in 1:ncats) {
      ## only look at named columns
      for (j in 1:nvars) {
        ## and keep only non-padding/zeroed names
        ssplit <- ifelse( i==1, "", "/")
        if (data.all$xlengths[i] >= j) {
          xvar.labels[j] <-
            paste( xvar.labels[j], ssplit, data.all$xnames[[i]][j], sep="")
          xvec[j,i] <- 1
        } else {
          xvar.labels[j] <- paste( xvar.labels[j], ssplit, "NA", sep="")
        }
      }
    }
    rownames(xvec) <- xvar.labels

    # parameter equality constraints
    if (!is.null(equality)) {
      eqlen <- length(equality);
      ieq <- 1;
      for (i in 1:eqlen) {
        ieq <- ieq + 1;
        xynames <- get.xynames(equality[[i]], data);
        nynames <- length(xynames$ynames);
        yidx <- match(xynames$ynames, data.all$ynames) ;
        if (any(is.na(yidx))) {
          if (print.level >= 0) {
            print("response used in equality constraints is not in the model");
            cat("responses in model:", data.all$ynames, "\n");
            cat("responses in equality:", xynames$ynames, "\n");
            print("terminating multinomRob")
          }
          return(NULL);
        }
        xidx <- list();
        xidxlen <- 0;
        for (j in 1:nynames) {
          xidx[[j]] <- match(xynames$xnames[[j]], data.all$xnames[[ yidx[j] ]]);
          if (any(is.na(xidx[[j]]))) {
            if (print.level >= 0) {
              print("regressor used in equality constraints is not in the model");
              cat("response:", xynames$ynames[j]);
              cat("regressors in model:", data.all$xnames[[ yidx[j] ]], "\n");
              cat("regressors in equality:", xynames$xnames, "\n");
              print("terminating multinomRob")
            }
            return(NULL);
          }
          xidxlen <- xidxlen + length(xidx[[j]]);
        }
        xyvmatrix <- matrix(0, 3, xidxlen);
        kk <- 0;
        for (j in 1:nynames) {
          for (k in 1:length(xidx[[j]])) {
            kk <- kk + 1;
            xyvmatrix[1,kk] <- xidx[[j]][k] ;
            xyvmatrix[2,kk] <- yidx[j] ;
            xyvmatrix[3,kk] <- xvec[xidx[[j]][k], yidx[j]] ;
          }
        }
        if (any(xyvmatrix[3,] != 1)) {  # need to consolidate equality constraints
          if (any(xyvmatrix[3,1] != xyvmatrix[3,])) {
            # not a subset of previous constraints
            xyvu <- unique(xyvmatrix[3,]);
            xyvu <- xyvu[xyvu>1];
            xyvmin <- min(xyvu);  # code to consolidate on
            for (k in 1:length(xyvu)) {
              xvec[xvec == xyvu[k]] <- xyvmin;
              xyvmatrix <- xyvmatrix[,xyvmatrix[3,] != xyvu[k]];
            }
            if (dim(xyvmatrix)[2] > 0) { # parameters not previously constrained
              for (k in 1:(dim(xyvmatrix)[2])) {
                xvec[xyvmatrix[1,k], xyvmatrix[2,k]] <- xyvmin;
              }
            }
          }
        }
        else { # consolidation not needed
          for (k in 1:(dim(xyvmatrix)[2])) {
            xvec[xyvmatrix[1,k], xyvmatrix[2,k]] <- ieq;
          }
        }
      }      
      if (print.level > 0) {
        cat("\nEquality constraints among parameters (after consolidation):\n");
        nidxvals <- length(idxvals <- sort(unique(xvec[xvec>1])));
        for (i in 1:nidxvals) {
          cat("Equality constrained set", i, "\n");
          for (j in 1:nvars)  for (k in 1:ncats) {
            if (xvec[j,k]==idxvals[i]) {
              cat("outcome", data.all$ynames[k], "regressor",
                  data.all$xnames[[k]][j], "\n");
            }
          }
        }
      }
    }
    
    if (print.level > 0) {
      cat("\nYour Model (xvec):\n")
      print(xvec)
    }

    TotalY <- apply( Y, 1, sum)
    tvec  <- xvec;
    nvars.unique <- sum(xvec == 1) + length(unique(xvec[xvec>1]));
#    nvars.total  <- nvars

    ## ------------------- END: PARSING / MODEL BUILDING AUTOMATION ------------------

    #let's create jacstack
    jacstack  <- jacstack.function(X,nvars.unique,xvec)

    # check for regressors with a distinct value at only one observation
    jsingle <- jacstack.singles(jacstack);
    if (any(jsingle) & print.level > 0) {
      cat("\n");
      print("multinomRob:  WARNING.  Limited regressor variation...")
      print("WARNING.  ... A regressor has a distinct value for only one observation.")
      print("WARNING.  ... I'm using a modified estimation algorithm (i.e., preventing LQD")
      print("WARNING.  ... from modifying starting values for the affected parameters).")

      xsingle <- tvec != tvec;
      xsingle <- mnl.xvec.mapping(forward=FALSE,xvec,xsingle,jsingle,ncats,nvars);
      print("WARNING.  ... Affected parameters are TRUE in the following table.")
      cat("\n");
      print(xsingle);
      cat("\n");
    }

    #load up genoudParms
    if(is.null(genoud.parms))
      genoud.parms  <- list(Domains=NULL);
    genoud.parms  <- genoudParms(genoud.parms)

    # initialize with terrible fit values
    mnl.fit  <- 9999999999;
    multinomT.fit  <- 9999999999;
    starting.fit   <- 9999999999;
    multinomT.foo <- NULL
    if (is.null(starting.values)) {
      starting.values <- vector(mode="numeric", length=nvars.unique);      
      if (print.level > 0)
        cat("\n\nmultinomRob(): Grouped MNL Estimation\n");
      Yp  <- Y/TotalY
      Yp[,ncats]  <- 1-apply(as.matrix(Yp[,1:(ncats-1)]),1,sum);
      mnl1 <- multinomMLE(Y=Y, Ypos=Ypos, Xarray=X, xvec=xvec, jacstack=jacstack,
                          xvar.labels=xvar.labels, choice.labels=choice.labels,
                          MLEonly=MLEonly,
                          print.level=print.level);

      if(!MLEonly)
        {
          starting.values <- mnl1$coeffvec
          mnl.fit <-
            fit.multinomial.C.lqd2(starting.values,X,Y,Ypos,xvec,tvec,ncats,nvars,
                                   nvars.unique,obs,TotalY);
        }
      
      if (print.level > 0) {
        if(!MLEonly)
          cat("MNL LQD Fit:",mnl.fit,"\n");        
        cat("MNL Estimates:\n");
        print(mnl1$coefficients);
        cat("\n");
        cat("MNL SEs:\n");
        print(mnl1$se);        
        cat("\n");
      }

      if(MLEonly)
        {
          z <- list(coefficients=mnl1$coefficients,
                    se=mnl1$se,
                    mnl=mnl1,
                    fitted.values=mnl1$fitted.prob,
                    deviance=as.double(mnl1$GNlist$LLvals)*2,
                    value=as.double(mnl1$GNlist$LLvals),
                    error=mnl1$error,
                    type="MLEonly"
                    )
          class(z)  <- "multinomRob"          
          return(z)
        }

      starting.values.user  <- starting.values
      tvec    <- mnl.xvec.mapping(forward=FALSE,xvec,tvec,starting.values,
                                  ncats,nvars)
      starting.fit <-
         fit.multinomial.C.lqd2(starting.values,X,Y,Ypos,xvec,tvec,ncats,nvars,
                                nvars.unique,obs,TotalY)

      use.mnl  <- TRUE;
      if (multinom.t > 0 & all(Ypos)) {
        if (print.level > 0) {        
          cat("\n\nmultinomRob(): Calculating multinomial-t starting values.\n");
        }

        Sres.raw <- res.std(Y,TotalY,mnl1$fitted.prob);
        multinomT.startvalues  <- list()
        multinomT.startvalues$beta  <- starting.values;
        multinomT.startvalues$Omega <- var(Sres.raw);
        #multinomT.startvalues$Omega <- var(mnl1$res[,1:(ncats-1)])/(obs-nvars.unique);
        multinomT.startvalues$df    <- 10;
        multinomT.foo  <- multinomT(Yp=Yp, Xarray=X, xvec=xvec, jacstack=jacstack,
                                    start=multinomT.startvalues, nobsvec=TotalY,
                                    fixed.df=multinom.t.df);

        for (ii in 0:6) {
          if(is.null(multinomT.foo$se$beta[1]) & is.na(multinom.t.df)) {
            if (print.level > 2)
              cat("Fixing Bad Estimate:",ii,"\n")
            multinomT.foo  <- multinomT(Yp=Yp, Xarray=X, xvec=xvec, jacstack=jacstack,
                                        start=multinomT.startvalues, nobsvec=TotalY,
                                        fixed.df=(10^ii-.5));
            if (print.level > 2)
              cat("DF: ",multinomT.foo$par$df,"\n");    
            
            multinomT.startvalues$beta  <- multinomT.foo$par$beta
            multinomT.startvalues$Omega <- multinomT.foo$par$Omega
            multinomT.startvalues$df    <- multinomT.foo$par$df
            multinomT.foo  <- multinomT(Yp=Yp, Xarray=X, xvec=xvec, jacstack=jacstack,
                                       start=multinomT.startvalues, nobsvec=TotalY);
          } else if (is.null(multinomT.foo$se$beta[1]) & !is.na(multinom.t.df)) {
            if (print.level > 1)
              {
                cat("WARNING: Multinom-T SEs are null, but we are using multinomT point estimates for starting values anyways because you explicitly sepecifed a DF\n")
              }
            use.mnl  <- FALSE;
            starting.values <- multinomT.foo$par$beta;
            multinomT.fit <-
              fit.multinomial.C.lqd2(starting.values,X,Y,Ypos,xvec,tvec,ncats,nvars,
                                     nvars.unique,obs,TotalY);
            if (print.level > 1)
              {
                cat("Multinom-T LQD Fit (step ",ii,"):",multinomT.fit,"\n")  
                beta.print <-
                  mnl.xvec.mapping(forward=FALSE,xvec,tvec,multinomT.foo$par$beta,ncats,nvars)
                cat("Multinom-T Beta Estimates (step ",ii,"):\n")
                print( beta.print )
                cat("Multinom-T SEs are not defined\n")
              }
            break;
          } else  {
            use.mnl  <- FALSE;
            starting.values <- multinomT.foo$par$beta;
            multinomT.fit <-
              fit.multinomial.C.lqd2(starting.values,X,Y,Ypos,xvec,tvec,ncats,nvars,
                                     nvars.unique,obs,TotalY);
            if (print.level > 1)
              {
                cat("Multinom-T LQD Fit (step ",ii,"):",multinomT.fit,"\n")  
                beta.print <-
                  mnl.xvec.mapping(forward=FALSE,xvec,tvec,multinomT.foo$par$beta,ncats,nvars)
                se.print <-
                  mnl.xvec.mapping(forward=FALSE,xvec,tvec,multinomT.foo$se$beta,ncats,nvars)
                cat("Multinom-T Beta Estimates (step ",ii,"):\n")
                print( beta.print )
                cat("Multinom-T Beta SEs (step ",ii,"):\n")
                print( se.print )
                cat("Multinom-T Omega Estimates (step ",ii,"):\n")
                print( multinomT.foo$par$Omega )
                cat("Multinom-T DF (step ",ii,"):", multinomT.foo$par$df,"\n\n")
              }
            break;
          }
        }#ii
      }#end of multinom.t

      if (use.mnl==FALSE & !is.null(starting.values) & multinom.t < 2) {
        if (multinomT.fit >= mnl.fit)
          use.mnl <- TRUE;
      }

      if(use.mnl==TRUE) {
        if (print.level > 0)
          cat("multinomRob(): Using grouped MNL estimates as starting values.\n")
        
        starting.values <- vector(mode="numeric", length=nvars.unique);
        starting.values <-
          mnl.xvec.mapping(forward=TRUE,xvec,mnl1$coefficients,starting.values,
                           ncats,nvars);
      } else {
        if (print.level > 0) {
          cat("multinomRob(): Using multinomial-t estimates as starting values.\n")
          cat("Multinom-T LQD Fit:",multinomT.fit,"\n")  
          beta.print <-
            mnl.xvec.mapping(forward=FALSE,xvec,tvec,multinomT.foo$par$beta,ncats,nvars)
          cat("Multinom-T Estimates:\n")
          print( beta.print )
          cat("Multinom-T DF:", multinomT.foo$par$df,"\n\n")
        }
        starting.values <- multinomT.foo$par$beta        
      }
    }#is.null(starting.values)
      
    tvec <- mnl.xvec.mapping(forward=FALSE,xvec,tvec,starting.values,ncats,nvars);

    fit.starting <-
      fit.multinomial.C.lqd2(starting.values,X,Y,Ypos,xvec,tvec,ncats,nvars,
                             nvars.unique,obs,TotalY)

    if (print.level > 2) {
      cat("multinomRob(): Starting Values \n")
      print(tvec)
      cat("multinomRob(): starting fit =",fit.starting,"\n");
    }
        
    genoud.out  <- genoudRob(genoud.fun,
                          nvars=nvars.unique,
                          starting.values = starting.values,
                          as.list(genoud.parms))
      
    s0    <- genoud.out$value;
    lqd.beta.vector  <- genoud.out$par;
    tvec  <- xvec;
    tvec  <- mnl.xvec.mapping(forward=FALSE,xvec,tvec,lqd.beta.vector,ncats,nvars);

    if (print.level>0) {    
      cat("\n");
      cat("\nLQD Results:\n");
      print(tvec);
      cat("\nLQD sigma:",s0,"\n\n")
    }
    
    #############################################
    ## TANH                                     #
    #############################################
    
    cat("\n(multinomTanh):\n");
      
    # modify LQD results for single unique value regressors
    tvec  <- mnl.xvec.mapping(forward=FALSE,xvec,tvec,
                              ifelse(jsingle, starting.values, lqd.beta.vector),
                              ncats,nvars);

    mout <- multinomTanh(Y=Y, Ypos=Ypos, X=X, jacstack=jacstack,
                         xvec=xvec, tvec=as.data.frame(tvec), 
                         pop=TotalY,s2=s0^2,
                         xvar.labels=xvar.labels, choice.labels,
                         print.level=print.level);
    error <- mout$mtanh$error;        
    if (mout$mtanh$error==0) {
      if (print.level > 0) {
        cat("Tanh Estimates\n")
        print(mout$mtanh$coef)
        cat("\nTanh Sandwich SEs\n");
        print(mout$mtanh$se);
        cat("\n")
        cat("TANH sigma:",sqrt(mout$mtanh$tanhsigma2),"\n\n")        
      } #end print      
    } else {
      if (print.level >= 0) {
        cat("WARNING: Tanh returned an error code:",mout$mtanh$error,"\n");
        cat("WARNING: Optimization not complete.  Starting values are probably not good enough. \n")
        if(iter) {
          cat("WARNING: Hopefully we will fix this problem during subsequent iterations.\n")
        } else {
          cat("WARNING: Rerun with new different starting value or with the iteration option left on\n")
        }
      }
    }

    # If the fit value gets worse, we stop.  Which is *NOT*
    # what the reference file does.  The reference file
    # keeps going (even if the new value is worse until there is
    # stability.
    s.count  <- 0;
    if (iter)
      {
        s0.old  <- s0;
        s0.tanh.old  <- sqrt(mout$mtanh$tanhsigma2);
        tolerance  <- genoud.parms$solution.tolerance;

        #let's shrink our domains by a half for the iteration stuff because we are
        # going to assume that the tanh starting values are not that bad.
        genoud.parms$scale.domains  <- genoud.parms$scale.domains/2;

        mout.old  <- mout;
      }#iter
    while(iter)
      {
        s.count  <- s.count+1;
        if (print.level > 0) {
          cat("\n**********************************\n\n");
          cat("Iteration:",s.count,"\n")
        }

        if (!any(is.na(mout.old$mtanh$coeffvec))) {
          LQDfit.tanhcoefs  <-
            fit.multinomial.C.lqd2(mout.old$mtanh$coeffvec,X,Y,Ypos,xvec,tvec,
                                   ncats,nvars,nvars.unique,obs,TotalY);
        } else {
          LQDfit.tanhcoefs  <- 9999999999;
        }

        best  <- order(c(starting.fit, mnl.fit, multinomT.fit, LQDfit.tanhcoefs))[1]
        if (s.count > 1 && !any(is.na(mout.old$mtanh$coeffvec))) {
          if (print.level > 0)
            cat("using LQD sigma provided by Tanh.\n")
          starting.values  <- mout.old$mtanh$coeffvec
        } else if (best==1) {
          if (print.level > 0)
            cat("using best non-genoud LQD sigma. Provided by user starting values.\n")
          starting.values  <- starting.values.user
        } else if (best==2) {
          if (print.level > 0)
            cat("using best non-genoud LQD sigma. Provided by ML multinomial.\n")
          starting.values  <- mnl1$coeffvec
        }  else if (best==3) {
          if (print.level > 0)
            cat("using best non-genoud LQD sigma. Provided by multinomial-t.\n")
          starting.values  <- multinomT.foo$par$beta;
        } else if (best==4) {
          if (print.level > 0)
            cat("using best non-genoud LQD sigma. Provided by Tanh.\n")
          starting.values  <- mout.old$mtanh$coeffvec
        }

        genoud.out  <- genoudRob(genoud.fun,
                                 nvars=nvars.unique,
                                 starting.values = starting.values,
                                 as.list(genoud.parms))
        
        s0   <- genoud.out$value;
        lqd.beta.vector  <- genoud.out$par;
        tvec  <- xvec;
        tvec  <- mnl.xvec.mapping(forward=FALSE,xvec,tvec,lqd.beta.vector,
                                  ncats,nvars);
        
        if (print.level>0) {    
          cat("\n(iter",s.count,"): LQD Results:\n");
          print(tvec);
          cat("\n(iter",s.count,"): LQD sigma:",s0,"\n\n")
        }

        if (s0 > (s0.old+tolerance))
            {
              if (print.level > 0)
                {
                  cat("WARNING: New fit is worse.  Stopping.\n");
                  cat("It may be a good idea to restart with a larger GENOUD population.\n");
                } #end of print.level
              error  <- -1;
              iter  <- FALSE;              
              break;
            }

        if (!is.na(s0.tanh.old) && converged.test(s0.old, s0, tolerance))
          {
              if (print.level > 0)
                {
                  cat("\n**********************************\n")
                  cat("CONVERGED\n");

                  if (print.level > 0) {
                    cat("(Converged): Tanh Estimates\n")
                    print(mout$mtanh$coef)
                    cat("\n(Converged): Tanh Sandwich SEs\n");
                    print(mout$mtanh$se);
                    cat("\n")
                    cat("(Converged): TANH sigma:",sqrt(mout$mtanh$tanhsigma2),"\n\n")
                    s0.tanh  <- sqrt(mout$mtanh$tanhsigma2);
                  } #end print                                    
                } #end of print.level
              error  <- 0;
              iter  <- FALSE;              
              break;            
          } else {
            #new s0 is smaller
        
            # modify LQD results for single unique value regressors
            tvec  <- mnl.xvec.mapping(forward=FALSE,xvec,tvec,
                                 ifelse(jsingle, starting.values, lqd.beta.vector),
                                 ncats,nvars);

            mout <- multinomTanh(Y=Y, Ypos=Ypos, X=X, jacstack=jacstack,
                                 xvec=xvec, tvec=as.data.frame(tvec), 
                                 pop=TotalY,s2=s0^2,
                                 xvar.labels=xvar.labels, choice.labels,
                                 print.level=print.level);
            
            error  <- mout$mtanh$error;        
            if (mout$mtanh$error==0)
              {
                if (print.level > 0) {
                  cat("(iter",s.count,"): Tanh Estimates\n")
                  print(mout$mtanh$coef)
                  cat("\n(iter",s.count,"): Tanh Sandwich SEs\n");
                  print(mout$mtanh$se);
                  cat("\n")
                  cat("(iter",s.count,"): TANH sigma:",sqrt(mout$mtanh$tanhsigma2),"\n\n")
                  s0.tanh  <- sqrt(mout$mtanh$tanhsigma2);

                  mout.old  <- mout;                  
                  s0.tanh.old  <- s0.tanh;
                } #end print                  
              } else {
                if (print.level > 0)
                  cat("(iter",s.count,"): tanh returns an error:",mout$mtanh$error,"\n");
              }
          } #end of else
            
        if (s.count==maxiter)
          {
            if (print.level > 0) 
              cat("\nWARNING: Maximum number of iterations reached.  Results have NOT converged\n\n");

            #iter failed
            error  <- -1;
            iter  <- FALSE;
            break;
          }
        if (s0 < s0.old)
          s0.old  <- s0;        
      }#iter check

    if (!exists("mnl1"))
      mnl1  <- NULL;
    if (!exists("multinomT.foo"))
      multinomT.foo  <- NULL;

    #Rotated Residuals.  This will result in a relatively easy to interpret
    #vector of residuals.  But this residual vector will NOT be a
    #consistent set of ortho residuals.
    residuals.rotate  <- matrix(nrow=obs,ncol=ncats)
    if (mout$mtanh$error < 32) {
      for (ii in 1:ncats)
        {
          tindx  <- 1:ncats
          tindx[1]  <- ii;
          tindx[ii] <- 1;
          
          YTmp     <- Y[,c(tindx)]
          YposTmp  <- Ypos[,c(tindx)]
          XTmp     <- X[,,c(tindx), drop=FALSE];
          jacstackTmp  <- jacstack[,,c(tindx), drop=FALSE]
          tvec  <- as.matrix(mout$mtanh$coefficients[,c(tindx), drop=FALSE])
          
          foo  <- permute(Y=YTmp, Ypos=YposTmp, Xarray=XTmp,
                          jacstack=jacstackTmp, tvec=tvec,
                          pop=TotalY, sigma=sqrt(mout$mtanh$disp),
                          weight=mout$mtanh$w)
          
          residuals.rotate[,ii]  <- foo$student[,1];
        }  #end of ii loop
    }
    
    z  <- list(coefficients=mout$mtanh$coefficients,
               se=mout$mtanh$se,
               LQDsigma2=mout$mtanh$dispersion,
               TANHsigma2=mout$mtanh$tanhsigma2,
               weights=mout$weights,
               Hdiag=mout$Hdiag,
               prob=mout$mtanh$prob,
               residuals.rotate=residuals.rotate,
               residuals.student=mout$cr$student,
               residuals.standard=mout$cr$standard,
               mnl=mnl1,
               multinomT  = multinomT.foo,
               genoud = genoud.out,
               mtanh = mout$mtanh,
               error = error,
               iter = s.count,#the iter number at the solution
               type="robust"); 
    class(z)  <- "multinomRob"
    return(z)
  }


#mimicking the definition of summary.default
#summary.default <-
#    function(object, ..., digits = max(3, getOption("digits") - 3))
summary.multinomRob <- function(object, ..., digits=3, weights=FALSE)
{
  
  
  if (class(object) != "multinomRob") {
    warning("Object not of class 'multinomRob'")
    return(NULL)
  } 

  if(object$type!="MLEonly")
    {
      out.mtanh <- list()
      
      if (!is.null(object$mtanh) && is.list(object$mtanh) && object$mtanh$error == 0) {
        ncats <- NCOL(object$mtanh$coefficients);
        nx <- NROW(object$mtanh$coefficients);
                                        # find length of longest regressor label
        labmax <- 0;
        for (i in 1:ncats) {
          for (j in 1:nx) {
            labn <- nchar(strsplit(row.names(object$mtanh$se)[j], "/")[[1]][i]) ;
            if (labmax < labn) labmax <- labn;
          }
        }
                                        # generate a blank variable with as many spaces as the longest label
        spc <- " ";
        blank <- "";
    for (i in 1:labmax) blank <- paste(spc, blank, sep="");
        
        for (i in 1:ncats) {
          tmp <- as.data.frame(
                               list("Est" = object$mtanh$coefficients[,i],
                                    "SE Sand" = object$mtanh$se[,i],
                                        #                    "SE OPG"  = object$mtanh$se.opg[,i],
                                        #                    "SE Hess" = object$mtanh$se.hes[,i],
                                    "t-val Sand" = object$mtanh$coefficients[,i]/object$mtanh$se[,i]))
          tmp <- as.matrix(tmp);
          rn <- rep(blank, nx);
          for (j in 1:nx) {
            rnj <- strsplit(row.names(object$mtanh$se)[j], "/")[[1]][i] ;
            substr(rn[j], 1, nchar(rnj)) <- rnj;
          }
          dimnames(tmp)[[1]] <- rn;
          out.mtanh[[i]] <- tmp
          choice.lables  <- labels(object$coefficients)[[2]]
          cat("\nChoice",i,":",choice.lables[i],"Estimates and SE:\n")
          print(signif(tmp,digits=digits))
          cat("\n")
        }
        
        cat("\n");
        cat("LQD sigma:",sqrt(object$LQDsigma2),"\n")
        cat("TANH sigma:",sqrt(object$TANHsigma2),"\n\n")
        
        indx  <- object$weights[,2:ncats]==0;
        indx[is.na(indx)] <- FALSE;
        if (ncats==2) indx <- matrix(indx, length(indx), 1);
        w0.obs    <- sum(apply(indx,1,any));
        cat("Number of Observations:",NROW(object$prob),"\n")
        cat("Number of observations with at least one zero weight:",w0.obs,"\n")
        w0  <- sum(indx)
        cat("Number of zero weights:",w0,"\n")
        
        if (weights) {
          cat("TANH: weights\n");
          print(object$weights);
        }
      }
      else {
        stop("error encountered in tanh result object");
      }
    } else { #END OF !MLEonly

      if (!is.null(object$mnl) && is.list(object$mnl) && object$mnl$error == 0) {
        ncats <- NCOL(object$mnl$coefficients);
        nx <- NROW(object$mnl$coefficients);
                                        # find length of longest regressor label
        labmax <- 0;
        for (i in 1:ncats) {
          for (j in 1:nx) {
            labn <- nchar(strsplit(row.names(object$mnl$se)[j], "/")[[1]][i]) ;
            if (labmax < labn) labmax <- labn;
          }
        }
                                        # generate a blank variable with as many spaces as the longest label
        spc <- " ";
        blank <- "";
    for (i in 1:labmax) blank <- paste(spc, blank, sep="");
        
        for (i in 1:ncats) {
          tmp <- as.data.frame(
                               list("Est" = object$mnl$coefficients[,i],
                                    "SE Sand" = object$mnl$se[,i],
                                        #                    "SE OPG"  = object$mnl$se.opg[,i],
                                        #                    "SE Hess" = object$mnl$se.hes[,i],
                                    "t-val Sand" = object$mnl$coefficients[,i]/object$mnl$se[,i]))
          tmp <- as.matrix(tmp);
          rn <- rep(blank, nx);
          for (j in 1:nx) {
            rnj <- strsplit(row.names(object$mnl$se)[j], "/")[[1]][i] ;
            substr(rn[j], 1, nchar(rnj)) <- rnj;
          }
          dimnames(tmp)[[1]] <- rn;
          choice.lables  <- labels(object$coefficients)[[2]]
          cat("\nChoice",i,":",choice.lables[i],"Estimates and SE:\n")
          print(signif(tmp,digits=digits))
          cat("\n")
        }
        
        cat("\n");
        cat("Residual Deviance:",as.double(object$deviance),"\n\n")
      } else {
        stop("error encountered in MNL result object");      
      }
    }#end of MLEonly
} #end of summary.multinomRob()


#mimicking the definition of plot 
plot.multinomRob  <- function(x, ...)
  {
    #this plots the residuals
    if (class(x) != "multinomRob") {
      warning("Object not of class 'multinomRob'")
      return(NULL)
    } 

#    out.mtanh <- list()

    if (!is.null(x$mtanh) && is.list(x$mtanh) && x$mtanh$error == 0)
      {
        ncats <- NCOL(x$mtanh$coefficients);        
        choice.lables  <- labels(x$coefficients)[[2]]

        ask.start  <- par("ask",no.readonly=TRUE)

        for (i in 1:(ncats-1))
          {
            plot(x$residuals.student[,i],xlab="Observation", ylab="Residual", ...)
            title(main=paste(choice.lables[i],"\nStudentized Residuals"));

            if (i==1)
              par(ask=TRUE);
          } #end for
        par(ask=ask.start)
      }
    else {
      print("error encountered in tanh result object");
    }
    #end if
  }#end of plot.multinomRob


#Function to create ortho-residuals (with base correction) for the
#other choices.  This will result in a relatively easy to interpret
#vector of residuals.  But this residual vector will NOT be a
#consistent set of ortho residuals.
#
#Modeled on
#lapo:~/xchg/election/R/multinomial/FL2/base4.permute4.origtanh2.R
#permute.newtanh4.R
#
permute  <- function(Y, Ypos, Xarray, jacstack, tvec, pop, sigma, weight)
  {
#tvec: tanh estimated values    
    #from multinomTanh
    nobs  <- dim(Y)[1]
    ncats <- dim(Y)[2]
    Hdiag <- robustified.leverage(tvec, Y, Ypos, Xarray, pop, ifelse(weight >0,1,0),jacstack);
#    w.Hdiag <- as.data.frame(matrix(c(as.vector(1:nobs),
#                                     signif(weight), signif(Hdiag)),ncol=(ncats-1)+(ncats-1)+1));
#    names(w.Hdiag) <-
#      c("name",paste("weights:",choice.labels[1:ncats-1],sep=""),
#        paste("Hdiag:",choice.labels[1:ncats-1],sep=""));
    #cat("mtanh: weights, Hdiag (by choices)\n");

    cr <- fn.region.results(tvec, Y, Ypos, Xarray, pop, sigma, Hdiag);    
    return(list(pred=cr$pred, student=cr$student, standard=cr$standard, Hdiag=Hdiag))
  } #end of permute

.onAttach <- function( ... )
{
  MatchLib <- dirname(system.file(package = "multinomRob"))
  version <- packageDescription("multinomRob", lib.loc = MatchLib)$Version
  BuildDate <- packageDescription("multinomRob", lib.loc = MatchLib)$Date

  foo <- paste("## \n##  multinomRob (Version ", version, ", Build Date: ", BuildDate, ")\n",
               "##  See http://sekhon.berkeley.edu/robust for additional documentation.\n",
               "##  Please cite as: Walter R. Mebane, Jr. and Jasjeet S. Sekhon. \"Robust Estimation\n",
               "##   and Outlier Detection for Overdispersed Multinomial Models of Count Data\".\n",
               "##   American Journal of Political Science, 48 (April): 391-410. 2004.\n##\n",
               sep="")
  packageStartupMessage(foo)
}


