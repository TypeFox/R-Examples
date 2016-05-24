KSS.default <- function(formula,
                        additive.effects = c("none", "individual", "time", "twoways"),
                        consult.dim.crit = FALSE,
                        d.max            = NULL,
                        sig2.hat         = NULL,
                        factor.dim       = NULL,
                        level            = 0.01,
                        spar             = NULL,
                        CV               = FALSE,
                        convergence      = 1e-6,
                        restrict.mode    = c("restrict.factors","restrict.loadings"),
                        ...)
  {
    ##===================================================================================
    if(!class(formula)=="formula"){
      stop("\n Argument >>formula<< needs a formula-object like y~x1+... where the elements are matrices.")
    }
    if(!any(additive.effects==c("none", "individual", "time", "twoways"))){
      stop("\n Argument >>effect<< must be one of: \n none, individual, time, twoways")
    }
    if(!is.numeric(level)){
      stop("\n Argument >>alpha<< has to be numeric.")
    }
    if(!is.null(factor.dim) & !is.numeric(factor.dim)){
      stop("\n Argument >>factor.dim<< has to be numeric.")
    }
    if(!is.null(spar) & !is.numeric(spar)){
      stop("\n Argument >>spar<< has to be numeric.")
    }
    if(!is.logical(CV)){
      stop("\n Argument >>CV<< has to be TRUE or FALSE.")
    }
    if(!is.numeric(convergence)){
      stop("\n Argument >>convergence<< has to be numeric.")
    }
    if(!is.null(sig2.hat) &!is.numeric(sig2.hat)){
      stop("\n Argument >>sig2.hat<< has to be numeric.")
    }
    if(!any(restrict.mode==c("restrict.factors","restrict.loadings"))){
      stop("\n Argument >>restrict.mode<< must be either: restrict.factors or: restrict.loadings")
    }
    ##====================================================================================
    
    ## check "effect" and "dim.crit"
    effect        <- match.arg(additive.effects)
        
    ## extract data from formula
    names  <- names(model.frame(formula))
    PF.obj <- FUN.Pformula(formula = formula, effect = effect)
    
    N <- ncol(PF.obj[[1]]$ODM)
    T <- nrow(PF.obj[[1]]$ODM)
    P <- length(PF.obj)-1
    dat.dim 	  <- c(T, N, P)
    is.intercept  <- PF.obj[[1]]$I

    ## *OR*iginal *dat*a
    ORdat    <- sapply(1:(P+1), function(i) PF.obj[[i]]$ODM)
    ## *TR*ansformed *dat*a
    TRdat         <- sapply(1:(P+1), function(i) PF.obj[[i]]$TDM)
    dat.matrix    <- sapply(1:(P+1), function(i) PF.obj[[i]]$TDV) 

    Or.Y     <- ORdat[, 1,       drop = FALSE]					# (TN x 1)
    Or.X     <- ORdat[, 2:(P+1), drop = FALSE]				        # (TN x P)

    TR.Y     <- TRdat[, 1,       drop = FALSE]					# (TN x 1)
    TR.X     <- TRdat[, 2:(P+1), drop = FALSE]				        # (TN x P)
    
    TR.Y.mat <- matrix(TR.Y, T,     N)					        # (T x N)
    TR.X.mat <- matrix(TR.X, T, (N*P))						# (T x NP)


    #####################################################################################################    
    ## GCV Procedure to determine the smoothing parameter
    iterateGCV <- FUN.iterate.GCV(TR.Y.mat, TR.X.mat, N, T, P)
    spar.GCV   <- iterateGCV$PCA$spar
    if(is.null(spar)){
      spar.low <- spar.GCV*0.75
    }else{
      spar.low  <- spar
    }
    if(CV){
      spar.interval.max <- spar.GCV
      spar.CV           <- KSS.CV(kappa.interv=c(.Machine$double.eps, spar.interval.max),
                                  Y=TR.Y, X=TR.X, N=N, T=T, P=P, spar.dim.fit=spar.low, tol=convergence)$minimum
      cat("\n CV-Optimization converged.\n")
      spar.low <- spar.CV
    }
    #####################################################################################################
    t.seq <- seq(0, 1, length.out=T)
    ## Determination of the 'norder' argument of function 'smooth.Pspline()'
    if(T>=5){
      norder <- 2
    }
    if(T==4){
      norder <- 1
    }
    if(T<4){
      stop("Time-Dimension T of the panel data has to be at least T=4.")
    }
    TR.Y.mat.smth  <- smooth.Pspline(x = t.seq, y = TR.Y.mat,      norder = norder, spar   = spar.low)$ysmth #(T x N)    
    TR.X.mat.smth  <- smooth.Pspline(x = t.seq, y = TR.X.mat,      norder = norder, spar   = spar.low)$ysmth #(T x NP)
    TR.X.mat.smth2 <- smooth.Pspline(x = t.seq, y = TR.X.mat.smth, norder = norder, spar   = spar.low)$ysmth #(T x NP)
    
    ## calculate beta coefficents

    TR.Y.smth        <- matrix(TR.Y.mat.smth,  nrow= (N*T), ncol = 1)	       # (TN x 1)
    TR.X.smth        <- matrix(TR.X.mat.smth,  nrow= (N*T), ncol = P)	       # (TN x P)
    TR.X.smth2       <- matrix(TR.X.mat.smth2, nrow= (N*T), ncol = P)	       # (TN x P)

    t.TR.X.TR.X      <- crossprod(TR.X)       			               # (PxP)
    t.TR.X.TR.X.smth <- crossprod(TR.X, TR.X.smth)		               # (PxP)
    t.TR.X.TR.X.smth2<- crossprod(TR.X, TR.X.smth2)		               # (PxP)
    
    t.TR.X.TR.Y      <- crossprod(TR.X, TR.Y)     		               # (Px1)
    t.TR.X.TR.Y.smth <- crossprod(TR.X, TR.Y.smth)   		               # (Px1)
    
    bloc1            <- t.TR.X.TR.X - t.TR.X.TR.X.smth     		       # (PxP)
    bloc2            <- t.TR.X.TR.Y - t.TR.X.TR.Y.smth     	               # (Px1)
    ## common-Slope.Coefficients:
    inv.bloc1 <- solve(bloc1)
    com.slops.0      <- inv.bloc1%*%bloc2				       # (Px1)
    ## calculate first step residuals and estimate dimension of factor-structure
    Residu.mat       <- matrix((TR.Y - (TR.X %*% com.slops.0)), T, N)

    ## functional pca
    fpca.fit.obj     <- fpca.fit(Residu.mat, spar=spar.low)

    ## Estimation of Dimension
    dim.criterion    <- c("PC1", "PC2", "PC3", "BIC3","IC1", "IC2", "IC3",
                          "IPC1","IPC2", "IPC3",
                          "ABC.IC1", "ABC.IC2", 
                          "KSS.C",
                          "ED",
                          "ER", "GR")
    Opt.dim.Output   <- as.matrix(sapply(dim.criterion, function(dim.criterion){
      EstDim(dim.criterion, Obj=Residu.mat, d.max=d.max, factor.dim=factor.dim,
             sig2.hat=sig2.hat, level=level)[2]}))

    ## User-Interface--------------------------------------------------------------------------- 
    Opt.dim.Output.BaiNg          <- c(as.numeric(Opt.dim.Output[1:7,1]))
    names(Opt.dim.Output.BaiNg)   <- c("PC1","PC2","PC3","BIC3","IC1","IC2","IC3")
    Opt.dim.Output.Bai            <- c(as.numeric(Opt.dim.Output[8:10,1]))
    names(Opt.dim.Output.Bai)     <- c("IPC1","IPC2","IPC3")
    Opt.dim.Output.Alessi         <- c(as.numeric(Opt.dim.Output[11:12,1]))
    names(Opt.dim.Output.Alessi)  <- c("ABC.IC1", "ABC.IC2")
    Opt.dim.Output.KSS            <- c(as.numeric(Opt.dim.Output[13,1]))
    names(Opt.dim.Output.KSS)     <- c(" KSS.C")
    Opt.dim.Output.Onatski        <- c(as.numeric(Opt.dim.Output[11,1]))
    names(Opt.dim.Output.Onatski) <- c(" ED")
    Opt.dim.Output.RH             <- c(as.numeric(Opt.dim.Output[12:13,1]))
    names(Opt.dim.Output.RH)      <- c(" ER","GR")
    if(is.null(factor.dim) && consult.dim.crit){
      cat("-----------------------------------------------------------\n")
      cat("Results of Dimension-Estimations");
      cat("\n\n-Bai and Ng (2002):\n");       print(Opt.dim.Output.BaiNg, quote = FALSE, na.print="");
      cat("\n\n-Bai (2004):\n");              print(Opt.dim.Output.Bai, quote = FALSE, na.print="");
      cat("\n\n-Alessi et al. (2010):\n");    print(Opt.dim.Output.Alessi, quote = FALSE, na.print="");
      cat("\n-Kneip et al. (2012):\n");       print(Opt.dim.Output.KSS, quote = FALSE, na.print="");
      cat("\n-Onatski (2009):\n");            print(Opt.dim.Output.Onatski, quote = FALSE, na.print="");
      cat("\n-Ahn and Horenstein (2013):\n"); print(Opt.dim.Output.RH, quote = FALSE, na.print="");
      otp <- as.numeric(c(Opt.dim.Output.BaiNg,Opt.dim.Output.KSS,Opt.dim.Output.Onatski,Opt.dim.Output.RH))
      cat("\n")
      cat("-----------------------------------------------------------\n")
      ## Auxiliary Fct's
      myASK <- function(){
        cat("Please, choose one of the proposed integers: ")
        readLines(con = stdin(), n = 1)
      }
      is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
        if(is.na(x)){
          result <- FALSE
        }else{
          if(is.numeric(x) & abs(x - round(x)) < tol){
            result <- min(as.numeric(Opt.dim.Output[,1]))<=x & x<= max(as.numeric(Opt.dim.Output[,1]))
          }else{
            result <- FALSE
          }
        }
        return(result)
      }
      used.dim      <- suppressWarnings(as.numeric(myASK()))
      ## Checking the Comand-Line input:
      while(!is.wholenumber(used.dim)){
        cat("Wrong input! You have to select a particular dimension by giving a numeric integer value between ", min(as.numeric(Opt.dim.Output[,1]))," and ", max(as.numeric(Opt.dim.Output[,1])),"! \n")
        used.dim      <- suppressWarnings(as.numeric(myASK()))
      }
      cat("Used dimension of unobs. factor structure is:", used.dim,"\n")
      cat("-----------------------------------------------------------\n")
    }
    ##------------------------------------------------------------------------------------------
    if(!is.null(factor.dim)){
        used.dim <- factor.dim
      }
    if(is.null(factor.dim) && !consult.dim.crit){
      used.dim <- c(as.numeric(Opt.dim.Output[13,1]))# KSS.C
    }
    ## now: 'used.dim' is specified 
    if(used.dim > 0){
      factors       <- fpca.fit.obj$factors[,  1:used.dim, drop= FALSE]
      loadings      <- fpca.fit.obj$loadings[, 1:used.dim, drop= FALSE]
      factor.stract <- tcrossprod(factors, loadings)
      ## Eventually for later extensions: logical argument two.step in order
      ## to allow for a conditional estimation of the slope-parameters given
      ## the estimated dimension d. 
      two.step <- FALSE
      if(two.step){
        ## re-estimate beta=========================================================
        NEW.TR.Y.mat  <- TR.Y.mat - factor.stract
        NEW.TR.Y      <- as.vector(NEW.TR.Y.mat)
        beta          <- qr.solve(TR.X, NEW.TR.Y)
        beta          <- matrix(beta, P,1)
        ##==========================================================================
      }else{
        beta          <- com.slops.0
      }
    }
    if(used.dim == 0){# no fact.-struct
      factors       <- NULL
      loadings      <- NULL
      factor.stract <- matrix(0,T,N)
      beta          <- qr.solve(TR.X, TR.Y) # OLS
    }
    ## Built up the return object *est*
    
    FUN.add.eff.obj <- FUN.add.eff(PF.obj        = PF.obj,
                                   fpca.fit.obj  = fpca.fit.obj,
                                   beta.hat      = beta)

    ## degrees.of.freedom =======================================================================
    degrees.of.freedom <- (T*N - (N+T)*used.dim - P - 
                           is.intercept -
                           N*(effect == "individual" | effect == "twoways") - 
                           T*(effect == "time"       | effect == "twoways"))
    
    ## re-estimation of sig2.hat (Paper KSS Section 3.4)==========================================
    YOVc            <- PF.obj[[1]]$TRm$OVc
    XOVc            <- sapply(2:(P+1), function(i)PF.obj[[i]]$TRm$OVc)
    Or.Y_minus_YOVc <- Or.Y - YOVc
    Or.X_minus_XOVc <- Or.X - matrix(rep(XOVc, each=(N*T)), N*T, P)
    
    sig2.hat        <- sum((TR.Y - TR.X %*% beta - matrix(factor.stract, N*T, 1))^2)/degrees.of.freedom
    
    
    
    ## estimation of Variances ==================================================================
    if(used.dim > 0){
      ## estimation of beta-variance beta.V
      if(!two.step){
        beta.V <- sig2.hat * solve(bloc1) %*% (t.TR.X.TR.X + t.TR.X.TR.X.smth2 - 2*t.TR.X.TR.X.smth) %*% solve(bloc1)
      }
      if(two.step){
        beta.V <- sig2.hat * solve(crossprod(TR.X))
      }
      ## estimation of Intercept-variance    
      if(is.intercept){
        Intercept.V <- (sig2.hat + matrix(colMeans(Or.X),1,P) %*% beta.V %*% t(matrix(colMeans(Or.X),1,P)))
      }else{Intercept.V <- NULL}
    }
    if(used.dim == 0){
      if(!is.intercept){
        beta.V      <- sig2.hat * solve(t.TR.X.TR.X)
        Intercept.V <- NULL
      }else{        
        beta.V       <- sig2.hat * solve(crossprod(TR.X))
        Intercept.V  <- (sig2.hat + matrix(colMeans(Or.X),1,P) %*% beta.V %*% t(matrix(colMeans(Or.X),1,P)))

      }
    }
    ## Fitted Values =============================================================================
    fitted.values      <- matrix(c(rep(FUN.add.eff.obj$mu, T*N) + rep(FUN.add.eff.obj$beta.0, N) +
                                   rep(FUN.add.eff.obj$tau, each=T) + Or.X %*% beta)
                                 , T, N) + factor.stract

    ## Return ====================================================================================
    est                    <- vector("list")
    est$dat.matrix         <- dat.matrix
    est$formula            <- formula
    est$dat.dim            <- dat.dim
    est$slope.para         <- beta
    est$beta.V             <- beta.V
    est$names              <- names         #Names of: dependent variable and regressors
    est$is.intercept       <- is.intercept  #Intercept: TRUE or FALSE
    est$additive.effects   <- effect        #Additive-Effect-Type
    est$Intercept          <- FUN.add.eff.obj$mu     # Intercept
    est$Intercept.V        <- Intercept.V
    est$Add.Ind.Eff        <- FUN.add.eff.obj$tau    # Add. indiv. Effects
    est$Add.Tim.Eff        <- FUN.add.eff.obj$beta.0 # Add. time Effects (beta.0-function)
    est$unob.factors       <- factors
    est$ind.loadings       <- loadings
    est$unob.fact.stru     <- factor.stract
    est$used.dim           <- used.dim
    est$optimal.dim        <- Opt.dim.Output
    est$fitted.values      <- fitted.values
    est$orig.Y             <- matrix(Or.Y, T, N)
    est$residuals          <- est$orig.Y - est$fitted.values
    est$sig2.hat           <- sig2.hat
    est$degrees.of.freedom <- degrees.of.freedom
    est$call               <- match.call()
    class(est)             <- "KSS" 
    ##=============================================================================================
    return(est)
  }
