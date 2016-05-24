"tune.ahazpen"<-function(surv, X, weights, standardize = TRUE, penalty=lasso.control(), tune=cv.control(),dfmax=nvars,lambda,...)
  {
    ## Purpose: Tuning parameter selection for additive hazards lasso
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   surv       : Surv object (right censored/counting process)
    ##   X          : Numeric matrix of covariates
    ##   weights    : Weight vector (nonnegative)
    ##   standardize: Standardize X?
    ##   K          : Number of folds
    ##   trace      : Print out progress?
    ##   allfolds   : Optional user-specified list of folds
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen

    this.call <- match.call()

    nvars<-ncol(X)
    nobs<-nrow(X)
    
    if(!missing(lambda))
      if(length(lambda) < 2)
        stop("'lambda' should have length > 1")

    tune<-eval(tune)
    # Get tuning controls
    if(is.character(tune)){
      tmp<-c("cv.control","bic.control",tune)[pmatch(tolower(tune),c("cv.control","bic.control"),nomatch=3)]
      tune<-get(tmp,mode="function")
    }
    if(is.function(tune)) tune <- tune()
    if(is.null(tune$type)) {
      print(tune)
      stop("'tune' not recognized")
    }
    X<-as.matrix(X)
    
    # Check for weights
    if(missing(weights))
      weights<-rep(1,nrow(X))


    # Preliminary 'ahazpen' fit
    if(missing(lambda))
      fit <-ahazpen(surv=surv,X=X,weights=weights,standardize=standardize,penalty=penalty,dfmax=dfmax,...)
    else
       fit <-ahazpen(surv=surv,X=X,weights=weights,standardize=standardize,penalty=penalty,dfmax=dfmax,lambda=lambda,...)
    df <- fit$df

    # If no nonzero coefficients in fit, no point in continuing
    if(max(df)==0)
      stop("no nonzero coefficients")

    lambdanew <- fit$lambda
    nlambda <- length(lambdanew)
            
    # CROSS-VALIDATION
    if(tune$type=="CV")
      {
        foldsused<-list()
        k<-1
        for(j in 1:tune$rep)
          {
            if(tune$trace)
              cat(paste("Repetition: ",j,"/",tune$rep,"\n",sep=""))
            cvf<-tune$getfolds(nrow(X))
            foldsused[[k]]<-cvf;k<-k+1
            if(j==1){
              error <- matrix(0, nrow = cvf$nfolds, ncol = length(lambdanew))
              ll<-matrix(0,nrow=cvf$nfolds,ncol=tune$rep)
            }
            
            for(i in 1:cvf$nfolds)
              {
                if(tune$trace)
                  cat(paste("  Fold: ",i,"/",cvf$nfolds,"\n",sep=""))
                omit <- cvf$all.folds[[i]]

                # Turning warnings into errors - otherwise we cannot be sure that CV makes sense
                tryCatch(  tmp <- ahazpen(surv=surv[-omit,], X=X[-omit,],weights=weights[-omit],lambda=lambdanew,penalty=penalty,...),
                         warning = function(w) {
                           stop(paste("Caught fatal warning from ahazpen:",conditionMessage( w )),call.=FALSE)
                         })

                include <- ahaz.nzcoef(coef(tmp))
                beta<-as.matrix(coef(tmp)[include,,drop=FALSE])
                if (sum(include)>0) {
                  test <- ahaz(surv[omit,], X[omit,include],weights=weights[omit])
                  error[i,] <-error[i,]+apply(beta,2,function(x){ahaz.mse(test, x)})
                }
              }
          }
        tunem  <- apply(error/tune$rep,2,mean)
        tunesd <- apply(error/tune$rep,2,sd)
        lambda.min <- max(lambdanew[tunem<=min(tunem)])
        
        out <- list("lambda" = lambdanew, "tunem" = tunem,
                    "tunesd" = tunesd, "tuneup" = tunem + tunesd,
                    "tunelo" = tunem - tunesd, "lambda.min" = lambda.min,
                    "df" = df, "call" = this.call,
                    "tune"=tune,"penalty"=penalty,"foldsused"=foldsused,
                    "nfolds"=cvf$nfolds,"fit"=fit)
        
        class(out) <- "tune.ahazpen"
        return(out)
      } else if(tune$type=="BIC")
        {
          include <- ahaz.nzcoef(coef(fit))
          beta<-as.matrix(coef(fit)[include,,drop=FALSE])
          m<-ahaz(surv,X[,include],weights)
          factor<-tune$factor(sum(weights))
          invB<-ahaz.ginv(m$B)
          invD<-ahaz.ginv(m$D)
          kappa<-drop(t(m$d)%*%invB%*%m$d/(t(m$d)%*%invD%*%(m$d)))
          bic<-kappa * apply(beta,2,function(x){ahaz.mse(m,x)}) + fit$df * factor
          lambda.min <- max(lambdanew[bic<=min(bic)])
          out<-list("df"=fit$df,"lambda"=lambdanew,"tunem"=bic, "lambda.min"=lambdanew[which.min(bic)],
                     "call" = this.call,"tune"=tune,"penalty"=penalty,"fit"=fit)
          class(out) <- "tune.ahazpen"
          return(out)
        }  else {
          stop("invalid 'tune'")
        }
  }
