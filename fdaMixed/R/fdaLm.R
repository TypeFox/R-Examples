fdaLm <- function(formula,data,design,
                  boxcox=NULL,G=1,lambda=1,nlSearch=TRUE,
                  K.order=1,D.order=NULL,Fleft="tied",Fright="tied",
                  left=NULL,right=NULL
#                  ,lag.max=25
                  ){
  # --------------------------------------
  # Take call
  # --------------------------------------
  
  mf.call <- match.call(expand.dots=FALSE)

  # --------------------------------------
  # Extract observation and identifier
  # --------------------------------------
  
  # Extract observation vector
  mf         <- mf.call[c(1,match(c("formula","data"),names(mf.call),0))]
  mf[[1]]    <- as.name("model.frame")
  mf$formula <- formula(Formula::Formula(formula),lhs=1,rhs=0)
  mf$drop.unused.levels <- TRUE
  Ymat       <- eval(mf,parent.frame())[,1]
  Ntotal     <- length(Ymat)
  
  # Extract sample identifier
  if (length(Formula::Formula(formula))[1]==1) {
    # Only one sample
    MM <- 1
  } else {
    # Look for sample identifier inside 'data'
    mf         <- mf.call[c(1,match(c("formula","data"),names(mf.call),0))]
    mf[[1]]    <- as.name("model.frame")
    mf$formula <- formula(Formula::Formula(formula),lhs=2,rhs=0)
    mf$drop.unused.levels <- TRUE
    id.data    <- try(as.factor(eval(mf,parent.frame())[,1]),silent=TRUE)
    if (class(id.data)!="try-error") {
      id.names.data <- levels(id.data)[unique(id.data)]
      MM.data       <- length(id.names.data)
    }
    
    # Look for sample identifier inside 'design'
    mf         <- mf.call[c(1,match(c("formula","design"),names(mf.call),0))]
    mf[[1]]    <- as.name("model.frame")
    mf$formula <- formula(Formula::Formula(formula),lhs=2,rhs=0)
    mf$drop.unused.levels <- TRUE
    names(mf)[match("design",names(mf),0)] <- "data"
    id.design  <- try(as.factor(eval(mf,parent.frame())[,1]),silent=TRUE)
    if (class(id.design)!="try-error") {
      id.names.design <- levels(id.design)[unique(id.design)]
      MM.design       <- length(id.names.design)
    }

    # Validate identifiers
    if (class(id.data)=="try-error") {
      if (class(id.design)=="try-error") {
        stop("Sample identifier not found")
      } else {
        MM       <- MM.design
        id.names <- id.names.design
      }
    } else {
      if (class(id.design)=="try-error") {
        MM       <- MM.data
        id.names <- id.names.data
      } else {
        if (MM.data!=MM.design) stop("Mismatch between length of sample identifiers inside 'data' and 'design'")
        if (any(is.na(match(id.names.data,id.names.design)))) stop("Mismatch between names of sample identifiers inside 'data' and 'design'")
        MM       <- MM.data
        id.names <- id.names.data
      }
    }
  }

  # ------------------------------------
  # Validate observation and identifier
  # ------------------------------------
  
  NN <- Ntotal %/% MM  
  if ((Ntotal %% MM)!=0) stop("Observation length not a multiple of sample numbers")
  if (NN < 3) stop("Sample length must be at least 3")

  # Copy orignal if needed in Box-Cox analysis
  if (is.numeric(boxcox)) {
    if (!all(Ymat > 0)) stop("In Box-Cox analysis observations must be positive")
    geometricMean <- exp(mean(log(Ymat)))
    Yoriginal     <- matrix(Ymat,NN,MM)   # matrix() to make actual copy
  }

  # Reshape observation vector as a matrix
  attr(Ymat,"dim") <- c(NN,MM)

  # ------------------------------------------
  # Make design matrices for fixed effect
  # ------------------------------------------

  # Checks whether all variables are inside 'design' data frame
  # Might be invoked later!!!
  #if (is.element("design",names(mf.call))) {
  #  mf           <- mf.call[c(1,match("design",names(mf.call)))]
  #  mf[[1]]      <- as.name("colnames")
  #  names(mf)[2] <- "x"
  #  if (all(is.element(as.character(attr(terms(formula(Formula(formula),lhs=0,rhs=1)),"variables"))[-1],eval(mf,parent.frame())))) warning("All fixed effect design variables found in 'design' data frame")
  #}

  # Look for design variables inside 'design'
  look.in.data <- FALSE
  if (is.element("design",names(mf.call))) {
    mf         <- mf.call[c(1,match(c("formula","design"),names(mf.call),0))]
    mf[[1]]    <- as.name("model.frame")
    mf$formula <- formula(Formula::Formula(formula),lhs=0,rhs=1)
    mf$drop.unused.levels <- TRUE
    names(mf)[match("design",names(mf),0)] <- "data"
    mf         <- try(suppressWarnings(eval(mf,parent.frame())),silent=TRUE)
    if (class(mf)=="try-error") {
      warning("All fixed effect design variables not found in 'design' and/or parent environment. I will look in 'data' frame")
      look.in.data <- TRUE
    }
  } else {
    look.in.data <- TRUE
  }
  
  if (look.in.data) {
    # Not all design variables found inside 'design'
    # Look for design variables inside 'data'
    mf         <- mf.call[c(1,match(c("formula","data"),names(mf.call),0))]
    mf[[1]]    <- as.name("model.frame")
    mf$formula <- formula(Formula::Formula(formula),lhs=0,rhs=1)
    mf$drop.unused.levels <- TRUE
    mf         <- try(suppressWarnings(eval(mf,parent.frame())),silent=TRUE)
    if (class(mf)=="try-error") stop("All fixed effect design variables not found in 'data' and/or parent environment")
  }

  # If identifiers in both 'data' and 'design' with same length as data frame,
  # then reorder data frame  
  if ((MM > 1) && (class(id.data)!="try-error") && (class(id.design)!="try-error") && (dim(mf)[1]==MM)) {
    mf <- mf[match(id.names.data,id.names.design),,drop=FALSE]
  }

  # Make design matrix
  GammaMat             <- model.matrix(Formula::Formula(formula),data=mf,lhs=0,rhs=1)
  beta.names           <- colnames(GammaMat)
  p0                   <- length(beta.names)
  Nfixed               <- dim(GammaMat)[1]
  attributes(GammaMat) <- NULL
  attr(GammaMat,"dim") <- c(Nfixed,p0)

  # ------------------------------------------
  # Make design matrices for random effect
  # ------------------------------------------
  
  # Design matrix for random effects 
  if (length(Formula::Formula(formula))[2]==1) {
    # No random effects
    q0   <- 0
    Zmat <- matrix(0,MM,q0)
  } else {
    # Look for design variables inside 'design'
    look.in.data <- FALSE
    if (is.element("design",names(mf.call))) {
      mf         <- mf.call[c(1,match(c("formula","design"),names(mf.call),0))]
      mf[[1]]    <- as.name("model.frame")
      mf$formula <- formula(Formula::Formula(formula),lhs=0,rhs=2)
      mf$drop.unused.levels <- TRUE
      names(mf)[match("design",names(mf),0)] <- "data"
      mf         <- try(suppressWarnings(eval(mf,parent.frame())),silent=TRUE)
      if (class(mf)=="try-error") {
        warning("All fixed random design variables not found in 'design' and/or parent environment. I will look in 'data' frame")
        look.in.data <- TRUE
      }
    } else {
      look.in.data <- TRUE
    }
    
    if (look.in.data) {
      # Not all design variables found inside 'design'
      # Look for design variables inside 'data'
      mf         <- mf.call[c(1,match(c("formula","data"),names(mf.call),0))]
      mf[[1]]    <- as.name("model.frame")
      mf$formula <- formula(Formula::Formula(formula),lhs=0,rhs=2)
      mf$drop.unused.levels <- TRUE
      mf         <- try(suppressWarnings(eval(mf,parent.frame())),silent=TRUE)
      if (class(mf)=="try-error") stop("All random effect design variables not found in 'data' and/or parent frame")
    }

    # If identifiers in both 'data' and 'design' with same length as data frame,
    # then reorder data frame
    if ((MM > 1) && (class(id.data)!="try-error") && (class(id.design)!="try-error") && (dim(mf)[1]==MM)) {
      mf <- mf[match(id.names.data,id.names.design),,drop=FALSE]
    }

    # Make design matrix
    Zmat             <- model.matrix(Formula::Formula(formula),data=mf,lhs=0,rhs=2)
    u.names          <- colnames(Zmat)
    q0               <- length(u.names)    
    Nrandom          <- dim(Zmat)[1]
    attributes(Zmat) <- NULL
    attr(Zmat,"dim") <- c(Nrandom,q0)
  }

  # ---------------------------------------------
  # Interpretation of intercepts if necessary
  # ---------------------------------------------
  
  if ((p0==1) && (beta.names=="(Intercept)")) {
    if ((q0==0) || ((q0==1) && (u.names=="(Intercept)"))) {
      # Interpret fixed effect as global intercept
      Nfixed   <- Ntotal
      GammaMat <- matrix(1,Nfixed,1)
    } else {
      # Interpret fixed effect as the random effect
      Nfixed   <- Nrandom
      GammaMat <- matrix(1,Nfixed,1)
    }
  }
  
  if ((q0==1) && (u.names=="(Intercept)")) {
    if ((p0==0) || ((p0==1) && (beta.names=="(Intercept)"))) {
      # Interpret random effect as global intercept
      Nrandom <- Ntotal
      Zmat    <- matrix(1,Nrandom,1)
    } else {
      # Interpret fixed effect as the random effect
      Nrandom <- Nfixed
      Zmat    <- matrix(1,Nrandom,1)
    }
  }
  
  # ----------------------------------------------
  # Validate dimensions and resolve ANOVA's 
  # ----------------------------------------------

  if ((p0>0) && (q0>0) && (Nfixed!=Nrandom)) stop("Mismatch between rows in design matrices for fixed and random effects")

  if ((p0>0) && (!is.element(Nfixed,c(0,MM,Ntotal))))  stop("Mismatch between sample identifier and design frame for fixed effect")
  
  if ((q0>0) && (!is.element(Nrandom,c(0,MM,Ntotal)))) stop("Mismatch between sample identifier and design frame for random effect")
  
  anova.type <- "global"
  if ((p0>0) && (Nfixed==MM))  anova.type <- "marginal"
  if ((q0>0) && (Nrandom==MM)) anova.type <- "marginal"
    
  # ----------------------------
  # Interpret parameters
  # ----------------------------
  
  # Interpret lambda parameter
  lambda[lambda<0] <- 0
  coefsMat         <- matrix(0,sum(lambda>0),3)
  if (lambda[1]==0) stop("Leading lambda coefficient must be strictly positive")
  if (K.order %% 2 == 1) {
    if (length(lambda)>2) stop("For odd K.order lambda must have length 1 or 2")
    lambda.name   <- paste("-D^",2*K.order,sep="")
    coefsMat[1,3] <- -1
    if (length(lambda)>1) {
      lambda.name <- c(lambda.name,"I")
      if (lambda[2]>0) coefsMat[2,1] <- 1
    }
  } else {
    if (length(lambda) > 3) stop("For even K.order lambda must have length 1, 2 or 3")
    lambda.name   <- paste("D^",2*K.order,sep="")
    coefsMat[1,3] <- 1
    if (length(lambda)>1) {
      if (K.order %% 4 == 2) {
        lambda.name <- c(lambda.name,paste("-D^",K.order,sep=""))
        if (lambda[2]>0) coefsMat[2,2] <- -1
      } else {
        lambda.name <- c(lambda.name,paste("D^",K.order,sep=""))
        if (lambda[2]>0) coefsMat[2,2] <- 1
      }
    }
    if (length(lambda)>2) {
      lambda.name <- c(lambda.name,"I")
      if (lambda[3]>0) coefsMat[2+(lambda[2]>0),1] <- 1
    } 
  }

  # Initialize D.order is necessary
  if (is.null(D.order)) {
    D.order <- K.order
  } else {
    if ((nlSearch) && (D.order < K.order)) stop("Parameter search requested with D.order < K.order")
  }
  
  # Initialize boundary conditions
  if (is.character(Fleft) && (Fleft=="tied")) Fleft <- cbind(diag(K.order),matrix(0,K.order,K.order))
  if (is.character(Fleft) && (Fleft=="open")) Fleft <- cbind(matrix(0,K.order,1),diag(K.order),matrix(0,K.order,K.order-1))
  if (!is.matrix(Fleft)) stop("Fleft misspecified")
  if (any(dim(Fleft)!=c(K.order,2*K.order))) stop(paste("Matrix of left boundary conditions does not have dimension=(",K.order,",",2*K.order,")",sep=""))

  if (is.character(Fright) && (Fright=="tied")) Fright <- cbind(diag(K.order),matrix(0,K.order,K.order))
  if (is.character(Fright) && (Fright=="open")) Fright <- cbind(matrix(0,K.order,1),diag(K.order),matrix(0,K.order,K.order-1))
  if (!is.matrix(Fright)) stop("Fright misspecified")      
  if (any(dim(Fright)!=c(K.order,2*K.order))) stop(paste("Matrix of right boundary conditions does not have dimension=(",K.order,",",2*K.order,")",sep=""))

  # Initialize limits if necessary
  if (is.null(left))  left  <- 0
  if (is.null(right)) right <- NN
  if (left >= right) stop("Left limit must be strictly less than right limit")

  # ----------------------------
  # Initialize variables
  # ----------------------------
  
  # Initialize variables to contain results
  if (anova.type=="global") {
    attr(GammaMat,"dim") <- c(NN,MM*p0)
    attr(Zmat,"dim")     <- c(NN,MM*q0)
    Ymat                 <- cbind(Ymat,GammaMat,Zmat)
    projMat              <- matrix(0,NN,MM*(1+p0+q0)*(1+D.order))    
    betaHat              <- rep(0,p0)
    uBLUP                <- rep(0,q0)
  }
  if (anova.type=="marginal") {
    Yres              <- matrix(0,NN,MM)
    projMat           <- matrix(0,NN,MM*(1+D.order))
    betaHat           <- matrix(0,NN,p0)
    uBLUP             <- matrix(0,NN,q0)
    uBLUPinvG         <- matrix(0,NN,q0)
    proj.uBLUPinvG    <- matrix(0,NN,q0)
    condRes.uBLUPinvG <- rep(0,NN*q0)
    xBLUP.uBLUPinvG   <- matrix(0,NN*q0,1)
    logLik.uBLUPinvG  <- rep(0,5)
    dummyVec <- rep(0,0)
    dummyMat <- matrix(0,0,0)  
  }
  condRes   <- rep(0,NN*MM)
  xBLUP     <- matrix(0,NN*MM,1+D.order)
  logLikVec <- rep(0,5+D.order)
  Cbeta     <- matrix(0,p0,p0)
  Cu        <- matrix(0,q0,q0)
  
  # ---------------------------
  # Likelihood function
  # ---------------------------

  minus2logLik <- function(param,return.sigma=FALSE) {
    # -------------------------------------------------------------------------
    # Extract parameters:
    # Present implementation of variance parameters, might be improved later:
    #   1) If q0!=0 then G is assumed to be proportional to identity matrix
    #   2) Parameters modeled on log scale
    # -------------------------------------------------------------------------

    if (is.numeric(boxcox)) {
      mu <- param[1]
      if (q0 > 0) {
        G     <- exp(param[2])*diag(q0)
        coefs <- as.vector(exp(param[-c(1,2)])%*%coefsMat)
      } else {
        G     <- diag(0)
        coefs <- as.vector(exp(param[-1])%*%coefsMat)
      }
    } else {
      if (q0 > 0) {
        G     <- exp(param[1])*diag(q0)
        coefs <- as.vector(exp(param[-1])%*%coefsMat)
      } else {
        G     <- diag(0)
        coefs <- as.vector(exp(param)%*%coefsMat)
      }
    }

    # -----------------------------
    # Box-Cox transformation
    # -----------------------------
    
    if (is.numeric(boxcox)) {
      # Do Box-Cox transformation on untransformed observations
      .Call("boxcoxTransform",
            mu,
            geometricMean,
            Yoriginal,
            Ymat,
            DUP=FALSE,package="fdaMixed")
    }

    # --------------------
    # Find eigenvalues
    # --------------------
    
    roots <- findRoots(c(1,0,0)+(right-left)/NN*coefs,K.order)

    # ---------------------
    # Trace integral
    # ---------------------

    innerIntegral <- function(u) {
      # Find eigenvalues
      roots <- findRoots(cbind(u+(right-left)/NN*coefs[1],
                               (right-left)/NN*coefs[2],
                               (right-left)/NN*coefs[3]),K.order)

      # Make call
      return(as.vector(.Call("fdaTrace",
                             left,right,coefs[3],
                             Re(roots$left),Im(roots$left),
                             Re(roots$right),Im(roots$right),
                             Fleft,Fright,NN,
                             DUP=FALSE,PACKAGE="fdaMixed")))
    }
    
    # --------------------------------------------
    # Marginal ANOVA (multivariate calibration)
    # --------------------------------------------

    if (anova.type=="marginal") {
      # Call marginal ANOVA's
      .Call("marginalANOVA",
            G,
            Ymat,
            GammaMat,Zmat,
            Yres,
            betaHat,uBLUP,uBLUPinvG,
            Cbeta,Cu,
            DUP=FALSE,package="fdaMixed")
    
      # Compute squared length of projected random effects
      if (q0==0) {
        quad.uBLUP <- 0
      } else {
        # Call FDA projection
        .Call("fdaEngine",
              left,right,(right-left)/NN*coefs[3],
              Re(roots$left),Im(roots$left),Re(roots$right),Im(roots$right),
              Fleft,Fright,
              G,
              uBLUPinvG,
              proj.uBLUPinvG,dummyVec,dummyVec,xBLUP.uBLUPinvG,condRes.uBLUPinvG,
              dummyMat,dummyMat,logLik.uBLUPinvG,
              DUP=FALSE,PACKAGE="fdaMixed")

        # Second quadratic term
        quad.uBLUP <- sum(as.vector(uBLUP)*condRes.uBLUPinvG)
      }
    
      # ----------------------------------------------
      # Fit serial correlated effects
      # ----------------------------------------------

      # Make call
      .Call("fdaEngine",
            left,right,(right-left)/NN*coefs[3],
            Re(roots$left),Im(roots$left),Re(roots$right),Im(roots$right),
            Fleft,Fright,
            G,
            Yres,
            projMat,dummyVec,dummyVec,xBLUP,condRes,dummyMat,dummyMat,logLikVec,
            DUP=FALSE,PACKAGE="fdaMixed")

      # ---------------------------
      # Estimate noise variance
      # ---------------------------
      
      if (K.order %% 2 == 1) {
        sigma2hat <- sum(c(logLikVec[3],quad.uBLUP,
                            coefs[1]*logLikVec[5],
                           -coefs[3]*logLikVec[5+K.order])
                         )/(NN*(MM-p0))
      } else {
        sigma2hat <- sum(c(logLikVec[3],quad.uBLUP,
                            coefs[1]*logLikVec[5],
                           -coefs[2]*logLikVec[5+K.order/2],
                            coefs[3]*logLikVec[5+K.order])
                         )/(NN*(MM-p0))
      }

      # ---------------------------
      # Compute -2*log(likelihood)
      # ---------------------------

      res <- NN*(MM-p0)*log(sigma2hat)+
             (MM-p0)*integrate(innerIntegral,0,1,stop.on.error=FALSE)$value+
             NN*(MM-p0)
      if (q0!=0) res <- res + NN*determinant(diag(q0)+t(Zmat)%*%Zmat%*%G)$modulus
      if (p0!=0) res <- res - NN*determinant(Cbeta)$modulus
      attributes(res) <- NULL      
    }

    # ----------------------------------------------
    # Global ANOVA: Fit all effects simultaneously
    # ----------------------------------------------

    if (anova.type=="global") {
      # Make call
      .Call("fdaEngine",
            left,right,(right-left)/NN*coefs[3],
            Re(roots$left),Im(roots$left),Re(roots$right),Im(roots$right),
            Fleft,Fright,
            G,
            Ymat,
            projMat,betaHat,uBLUP,xBLUP,condRes,Cbeta,Cu,logLikVec,
            DUP=FALSE,PACKAGE="fdaMixed")

      # ---------------------------
      # Estimate noise variance
      # ---------------------------
      
      if (K.order %% 2 == 1) {
        sigma2hat <- sum(c(logLikVec[3:4],coefs[1]*logLikVec[5],
                           -coefs[3]*logLikVec[5+K.order])
                         )/(NN*MM-p0)
      } else {
        sigma2hat <- sum(c(logLikVec[3:4],coefs[1]*logLikVec[5],
                           -coefs[2]*logLikVec[5+K.order/2],
                            coefs[3]*logLikVec[5+K.order])
                         )/(NN*MM-p0)
      }
    

      # ---------------------------
      # Compute -2*log(likelihood)
      # ---------------------------

      res <- (NN*MM-p0)*log(sigma2hat)+
             MM*integrate(innerIntegral,0,1,stop.on.error=FALSE)$value+
             logLikVec[1]+logLikVec[2]+NN*MM-p0
    }

    # -------------------
    # Return result
    # -------------------
    
    if (return.sigma) {
      return(list(logLik=res,sigma2hat=sigma2hat))
    } else {
      return(res)
    }

    # End objective function: -2*log(likelihood)
  }

  # ----------------------------------------
  # Do non-linear optimization if requested
  # ----------------------------------------

  if (q0 > 0) {
    if (is.numeric(boxcox)) {
      param    <- c(boxcox[1],log(G),log(lambda[lambda>0]))
      if (nlSearch) {
        res    <- nlminb(param,minus2logLik)
        counts <- 1+sum(c(1,length(param))*res$evaluations)
        param  <- res$par
      } else {
        counts <- 1
      }
      boxcox           <- param[1]
      G                <- exp(param[2])*diag(q0)
      lambda[lambda>0] <- exp(param[-c(1,2)])      
    } else {
      param    <- c(log(G),log(lambda[lambda>0]))
      if (nlSearch) {
        res    <- nlminb(param,minus2logLik)
        counts <- 1+sum(c(1,length(param))*res$evaluations)
        param  <- res$par
      } else {
        counts <- 1
      }
      boxcox           <- "not done"
      G                <- exp(param[1])*diag(q0)
      lambda[lambda>0] <- exp(param[-1])
    }
  } else {
    if (is.numeric(boxcox)) {
      param    <- c(boxcox[1],log(lambda[lambda>0]))
      if (nlSearch) {
        res    <- nlminb(param,minus2logLik)
        counts <- 1+sum(c(1,length(param))*res$evaluations)
        param  <- res$par
      } else {
        counts <- 1
      }
      boxcox           <- param[1]
      G                <- diag(q0)
      lambda[lambda>0] <- exp(param[-1])
    } else {
      param    <- log(lambda[lambda>0])
      if (nlSearch) {
        res    <- nlminb(param,minus2logLik)
        counts <- 1+sum(c(1,length(param))*res$evaluations)
        param  <- res$par
      } else {
        counts <- 1
      }
      boxcox           <- "not done"
      G                <- diag(q0)
      lambda[lambda>0] <- exp(param)
    }
  }

  # Extract L-coefficients and estimate error variance
  # NB: Variable 'res' is used in final return, and hence should not be distorted
  res <- minus2logLik(param,return.sigma=TRUE)

  # -----------------------------------------
  # Compute conditional variance
  #  (experimental code, might be improved)
  # Non-correct formula used below. Hence excluded in Version 0.1. To be updated later.
  # -----------------------------------------

#  # Make coefficients
#  coefs <- lambda[lambda>0]%*%coefsMat
#
#  # Find variance of conditional residuals
#  if (anova.type=="global") {
#    # Diagonal of Green function corresponding to R*solve(I+R)
#    invAdiag <- rep(0,NN)
#    acfVec   <- rep(0,min(1+lag.max,NN-2))
#    roots <- findRoots(c(1,0,0)+(right-left)/NN*coefs,K.order)   
#    .Call("diagGreen",
#          left,right,(right-left)/NN*coefs[3],
#          Re(roots$left),Im(roots$left),Re(roots$right),Im(roots$right),
#          Fleft,Fright,
#          invAdiag,acfVec,
#          DUP=FALSE,PACKAGE="fdaMixed")
#
#    # Variance of conditional residuals
#    condResVar <- outer(invAdiag,rep(res$sigma2hat,MM))
#    acfVec <- acfVec/acfVec[1]
#    
#    warning("Variance of conditional residuals not fully implemented for global ANOVA")
#    warning("Auto covariance function might be wrong for global ANOVA")
#  }
#
#  if (anova.type=="marginal") {
#    # Diagonal of Green function corresponding to R*solve(I+R)
#    invAdiag <- rep(0,NN)
#    acfVec   <- rep(0,min(1+lag.max,NN-2))
#    roots <- findRoots(c(1,0,0)+(right-left)/NN*coefs,K.order)   
#    .Call("diagGreen",
#          left,right,(right-left)/NN*coefs[3],
#          Re(roots$left),Im(roots$left),Re(roots$right),Im(roots$right),
#          Fleft,Fright,
#          invAdiag,acfVec,
#          DUP=FALSE,PACKAGE="fdaMixed")
#
#    # Make matrix terms
#    if (q0==0) {
#      mat1 <- diag(MM)
#      mat2 <- diag(MM)
#    } else {
#      mat1 <- diag(MM)-Zmat%*%Cu%*%t(Zmat)
#      mat2 <- diag(MM)+Zmat%*%G%*%t(Zmat)
#    }
#    if (p0==0) {
#      mat3 <- diag(MM)
#    } else {
#      mat3 <- GammaMat%*%Cbeta%*%t(GammaMat)
#    }
#    
#    # Variance of conditional residuals
#    condResVar <- res$sigma2hat*acfVec[1]*
#      (mat1%*%mat2%*%mat1
#       -mat1%*%mat2%*%mat1%*%mat3%*%mat1
#       -mat1%*%mat3%*%mat1%*%mat2%*%mat1
#       +mat1%*%mat3%*%mat1%*%mat2%*%mat1%*%mat3%*%mat1)
#    if (MM > 0) {
#      colnames(condResVar) <- id.names
#      rownames(condResVar) <- id.names
#    }
#
#    # Rescale auto covariance as auto correlation
#    invAdiag <- invAdiag/acfVec[1]
#    acfVec   <- acfVec/acfVec[1]
#  }
  
  # --------------------
  # Reshape results
  # --------------------
  
  # Reshape xBLUP and conditional residuals
  attr(xBLUP,"dim")   <- c(NN,MM,1+D.order)
  attr(condRes,"dim") <- c(NN,MM)

  # Add appropriate name attributes
  names(lambda)          <- lambda.name
  if (p0 > 0) {
    if (anova.type=="global") {
      names(betaHat)     <- beta.names
    } else {
      colnames(betaHat)  <- beta.names
    }
    rownames(Cbeta)      <- beta.names
    colnames(Cbeta)      <- beta.names
  }
  if (q0 > 0) {
    if (anova.type=="global") {
      names(uBLUP)       <- u.names
    } else {
      colnames(uBLUP)    <- u.names
    }
    rownames(G)          <- u.names
    colnames(G)          <- u.names
    rownames(Cu)         <- u.names
    colnames(Cu)         <- u.names
  }
  if (MM > 1) {
    colnames(xBLUP)      <- id.names
    colnames(condRes)    <- id.names
#    colnames(condResVar) <- id.names                   # NB: Non-commented if condResVar is included
  }

  # --------------------
  # Return results
  # --------------------

#  if ((q0>0) && (p0>0)) warning("uVar is not corrected for estimation of fixed effect.\n  To be updated later.")
  
  return(list(logLik=res$logLik,ANOVA=anova.type,nlSearch=nlSearch,counts=counts,
              boxcoxHat=boxcox,Ghat=G,lambdaHat=lambda,sigma2hat=res$sigma2hat,
              betaHat=betaHat,uBLUP=uBLUP,xBLUP=xBLUP,
              condRes=condRes,
#              condResVar=condResVar,
#              invAdiag=invAdiag,acf=acfVec,
              betaVar=Cbeta
#              ,uVar=Cu
              ))
}
