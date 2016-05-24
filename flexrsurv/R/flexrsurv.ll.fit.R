flexrsurv.ll.fit<-function (X0, X, Z, Y, 
                            expected_rate, 
                            weights=NULL,
                            Spline_t0=BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                            Spline_t =BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t_NPH=TRUE,
                            bhlink=c("log", "identity"),
                            init=list(gamma0= NULL, alpha0=NULL, beta0=NULL, alpha=NULL, beta=NULL),
                            fastinit=TRUE,
                            optim.control=list(),
                            method=list(int_meth="CAV_SIM", step=diff(range(Y[,1]))/100, optim_meth="BFGS"),
                            vartype = c("oim", "opg", "none"), tol=1e-25,
                            debug=FALSE, debug.ll=FALSE, debug.gr=FALSE
                            )
{

# flexible relative survival model using full likelihood and 
# non iteratif, paramétrage identifiable
#
  #
  #         input :
  #
  #
  # gamma0 : vector of coef for baseline hazard
  # alpha0 ; vector of all coefs for non time dependant variables (may contain non-loglinear terms such as spline)
  # beta0 ; matrix of all coefs for log-linear but  time dependant variables  X%*%beta0(t)
  # beta  : matrix of coefs for beta(t) nTbasis * nTDvars for NLG and NPH
  # alpha : vector of coef for alpha(z) for NLG and NPH
  # X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
  # X : log lineair but time dependante variable 
  # Z : object of class time d?pendent variables (spline basis expended)
  # Y : object of class Surv
  # expected_rate : expected rate at event time T
  # weights : vector of weights  : LL = sum_i w_i ll_i
  #  Spline_t0, spline object for baseline hazard, with evaluate() méthod
  #  Intercept_t0=FALSE, option for evaluate, = TRUE all the basis, =FALSE all but first basis 
  #  Spline_t, spline object for time dependant effects,  with evaluate() méthod
  #  Intercept_t=FALSE, option for evaluate, = TRUE all the basis, =FALSE all but first basis 
  # init : list  of initial values
  # fastinit : if init=NULL, when fastinit=TRUE, init=(gamma0=rep(log(sum(status)/sum(time*(status==1)), ngamma0), othercoef=0)
  #                          when fastinit=FALSE, init in 3 steps: initgamma0, initalpha0alpah, initbeta0beta
  # optime.control : control parameters/options for optim()
  # method : optimisation method (optim_meth) for optim(), numerical intégration method (int_meth),
  # vartype : type of variance matrix : observed inf. mat (oim inv(-H)), robust/sandwich (robust H inv(S'S) H ),
  #           outer product of the gradients (opg inv(S'S)), wher where S is the matrix of scores

  # tol :   the tolerance for detecting linear dependencies in the columns of the information matrix
  # output : coef(with name, structure as attributes), var, conv & LL
  #

   bhlink  <- match.arg(bhlink)       # type baseline hazard

  vartype  <- match.arg(vartype)  # type of var
  dohessian <- ( vartype=="oim")


  time <- Y[, 1]
  status <- Y[, 2]

  n <- nrow(Y)

  # structure of parameter vector
  # baseline hazard 
   if(!is.null(Spline_t0)){
     nT0basis <- dim(evaluate(Spline_t0, time, intercept=Intercept_t0))[2]
     ngamma0 <- nT0basis
   }
   else {
     # no baseline hazard
     # the baseline is in NPH, nPHNLL or WCE effects
     nT0basis <- 0
     ngamma0 <- 0
   }
   
  # PHLIN and PHNLIN effects 
  if(!is.null(X0)){
    is.PH <- TRUE
    if(is.matrix(X0)) {
      nX0 <- dim(X0)[2]
    } else if(is.vector(X0)) {
      nX0 <- 1L
    } else {
      stop("wrong type of X0", call.=TRUE) 
    }
    nalpha0<-nX0
    ialpha0<-1:nX0 + ngamma0
    Ialpha0<-1:nX0 
    First.alpha0<-ngamma0+1
    first.alpha0<-1
  } else {
    is.PH <- FALSE
    nX0 <- 0L
    nalpha0 <- 0L
    alpha0 <- NULL
    ialpha0 <- NULL
    Ialpha0 <- NULL
    First.alpha0 <- NULL
    first.alpha0 <- NULL
  }
  
  # NPHLIN effects 
  if(!is.null(X)){
    is.NPHLIN <- TRUE
    nTbasis <- getNBases(Spline_t)
    nTbasis_NPH <- getNBases(Spline_t) - 1 + Intercept_t_NPH
    if(is.matrix(X)) {
      nX <- dim(X)[2]
    } else if(is.vector(X)) {
      nX <- 1L
    } else {
      stop("wrong type of X", call.=TRUE) 
    }
    nbeta0 <- sum(nTbasis_NPH)
    ibeta0 <- 1:nbeta0 + ngamma0 + nalpha0
    Ibeta0 <- 1:nbeta0 
    First.beta0 <- ngamma0 + nalpha0 + 1
    first.beta0 <- 1
  } else  {
    is.NPHLIN <- FALSE
    nTbasis <- 0L
    nTbasis_NPH <- 0L
    nX <- 0L
    nbeta0 <- 0L
    beta0 <- NULL
    ibeta0 <- NULL
    Ibeta0 <- NULL
    First.beta0 <- NULL
    first.beta0 <- NULL
  }

  
  # NPHNLIN effects 
  if(!is.null(Z)) {
    if( getNvar(Z)>0){
      is.NPHNLIN <- TRUE
      nTbasis <- getNBases(Spline_t)
      nZ<-getNvar(Z)

      nTbasis_NPHNLL <- rep(getNBases(Spline_t), nZ)
      nalpha <- getNparam(Z)
      ialpha <- 1:nalpha + ngamma0 + nalpha0 + nbeta0  
      Ialpha <- 1:nalpha  + nalpha0   
      First.alpha <- ngamma0 + nalpha0 + nbeta0 + 1  
      first.alpha <- nalpha0  + 1  
      
    # as first beta is constraints, nTbasis -1 beta per Z variable
      nbeta <- nZ * (nTbasis-1)
      ibeta <- 1:nbeta + ngamma0 + nalpha0 + nbeta0 + nalpha
      Ibeta <- 1:nbeta  + nbeta0 
      First.beta <- ngamma0 + nalpha0 + nbeta0 + nalpha + 1
      first.beta <- 1 + nbeta0 
    }
  } else  {
    is.NPHNLIN <- FALSE
#    nTbasis <- 0 already done with beta0
    nZ <- 0L
    nalpha <- 0L
    alpha <- NULL
    Ialpha <- NULL
    ialpha <- NULL
    First.alpha <- NULL
    first.alpha <- NULL
    nbeta <- 0L
    beta <- NULL
    Ibeta <- NULL
    ibeta <- NULL
    First.beta <- NULL
    first.beta <- NULL
  }

  # number of variables 
  nvar <- nX0 + nX + nZ
  # nb of parameters
  nparam = ngamma0 + nalpha0 + nbeta0 + nalpha + nbeta

  
  if (missing(expected_rate) || is.null(expected_rate)) 
    expected_rate <- rep(0, n)
  if (missing(weights) || is.null(weights)) {
    weights <- NULL
  } else {
    if (any(weights <= 0)) 
      stop("Invalid weights, must be >0", call.=TRUE)
    storage.mode(weights)    <- "double"
  }

# compute init values for gamma0
  lambda <- sum(status)/sum(time*(status==1))
  if( bhlink == "log"){
    gamma0init <- log(lambda)
  }
  else {
    gamma0init <- lambda
  }
  
  # control of init values
  if (!missing(init) && !is.null(init)) {
    if(ngamma0>0){
      if ( !is.null(init$gamma0)) {
        if (length(init$gamma0) != nT0basis){ 
          stop("Wrong length for initial values for gamma0", call.=TRUE)
        } else {
          initgamma0 <- init$gamma0
          do.init.gamma0 <- FALSE
        }
      } else {
        initgamma0 <- initcoef(Spline_t0, ncol=1L, init=gamma0init, intercept=Intercept_t0)
        do.init.gamma0 <- !fastinit
      }
    }
    else {
        initgamma0 <- NULL
    }
    if ( nX0 ) {
        if(!is.null(init$alpha0)) {
        if (length(init$alpha0) != nX0 ) {
          stop("Wrong length for initial values for alpha0", call.=TRUE)
        } else {
          initalpha0 <- init$alpha0
          do.init.alpha0 <- FALSE
          }
      } else {
          initalpha0 <- rep(0, nX0)
          do.init.alpha0 <- !fastinit
        }
      } else {
        initalpha0 <- NULL
        do.init.alpha0 <- FALSE
      }

    if (nX ){
        if (!is.null(init$beta0)) {
        if (length(init$beta0) != nbeta0){ 
          stop("Wrong length for initial values for beta0", call.=TRUE)
        } else {
          initbeta0 <- init$beta0
          do.init.beta0 <- FALSE
        }
      } else {
          initbeta0 <- rep(0, nbeta0)
          do.init.beta0 <-  !fastinit
        }
      } else {
      initbeta0 <- NULL
      do.init.beta0 <- FALSE
    }
    
    if (nZ ){
      if (!is.null(init$alpha)) {
        if (length(init$alpha) != nalpha ){ 
          stop("Wrong length for initial values for alpha", call.=TRUE)
        } else {
          initalpha <- init$alpha
          do.init.alpha <- FALSE
        }
      } else {
        initalpha <- rep(0, nalpha)
        do.init.alpha <-  !fastinit
      }        
      if (!is.null(init$beta)) {
        if (length(init$beta) != nbeta)  {
          stop("Wrong length for initial values for beta", call.=TRUE)
        } else {
          initbeta <- init$beta
          do.init.beta <- FALSE
        }
      } else {
        # if initalpha given (do.init.alpha=0), initbeta should be such that beta(t)=1
        # else initbeta=0
        initbeta <- rep(0, nbeta)
        do.init.beta <-  !fastinit
      }
    }
    else {
      initalpha <- NULL
      initbeta <- NULL
      do.init.alpha <- FALSE
      do.init.beta <- FALSE
    }        
# end if (!missing(init) && !is.null(init)) {
  } else {
    do.init <- !fastinit
    if(ngamma0>0){
      initgamma0 <- gamma0init * initcoef(Spline_t0, ncol=1L, intercept=Intercept_t0)
      optim.control.gamma0 <- optim.control
      optim.control.gamma0$parscale <- NULL
      optim.control.gamma0$ndeps <- NULL
      
      do.init.gamma0 <- !fastinit
    }else {
      initgamma0 <- NULL
      do.init.gamma0 <- FALSE
    }
    if ( nX0) {
      initalpha0 <- rep(0, nX0)
      do.init.alpha0 <- !fastinit
      } else {
      initalpha0 <- NULL
      do.init.alpha0 <- FALSE
    }
    if (nX) {
      initbeta0 <- rep(0, nbeta0)
      do.init.beta0 <- !fastinit
    } else {
      initbeta0 <- NULL
      do.init.beta0 <- FALSE
    }
    if (nZ) {
      initalpha <- rep(0, nalpha)
      initbeta <- rep(0 , nbeta)
      do.init.alpha <- !fastinit
      do.init.beta <- !fastinit
    } else {
      initalpha <- NULL
      initbeta <- NULL
      do.init.alpha <- FALSE
      do.init.beta <- FALSE
    }
  }# end else if (!missing(init) && !is.null(init)) {
  
    optim.control.gamma0 <- optim.control
    optim.control.gamma0$parscale <- NULL
    optim.control.gamma0$ndeps <- NULL
    optim.control.alpha <- optim.control
    optim.control.alpha$parscale <- NULL
    optim.control.alpha$ndeps <- NULL
    optim.control.beta <- optim.control
    optim.control.beta$parscale <- NULL
    optim.control.beta$ndeps <- NULL
  
  do.init <- do.init.gamma0 | do.init.alpha0 | do.init.beta0 | do.init.alpha | do.init.beta 
  
  storage.mode(initgamma0) <- "double"
  storage.mode(initalpha0) <- "double"
  storage.mode(initbeta0)  <- "double"
  storage.mode(initalpha)  <- "double"
  storage.mode(initbeta)   <- "double"

  
   
  # numerical integration method
  # computes steps for time integtration
  if(method$int_meth == "CAV_SIM"){
    int_meth <- "NC"
    intTD <- intTD_NC
    intTD_debug <- intTD_NC_debug
    intTD_base <- intTD_base_NC
    intTD_base_debug <- intTD_base_NC_debug
    intweightsfunc <- intweights_CAV_SIM
    step <- method$step
    mult <- 2
  } else if(method$int_meth == "SIM_3_8"){
    int_meth <- "NC"
    intTD <- intTD_NC
    intTD_base <- intTD_base_NC
    intweightsfunc <- intweights_SIM_3_8
    step <- method$step
    mult <- 3
  } else if(method$int_meth == "BOOLE"){
    int_meth <- "NC"
    intTD <- intTD_NC
    intTD_base <- intTD_base_NC
    intweightsfunc <- intweights_BOOLE
    step <- method$step
    mult <- 4      
  } else if(method$int_meth == "GLM"){
    int_meth <- "GLM"
    intTD <- intTD_GLM
    intTD_base <- fastintTD_base_GLM
    intweightsfunc <- NULL
    step <- GLMStepParam(cuts=method$bands)
    Nstep <- WhichBand(time, step)-1L
  }

  if( int_meth == "NC"){
    STEPS <- cutT(time, step=method$step, mult=mult)
    Nstep <- STEPS$NstepT
    step <- STEPS$stepT
  }

   
    LL <- +Inf

# objective, gradiant functions
  if( bhlink == "log"){
    ll_G0A0B0AB    <- ll_flexrsurv_GA0B0AB
    ll_gamma0      <- ll_flexrsurv_gamma0
    ll_alpha0alpha <- ll_flexrsurv_alpha0alpha
    ll_beta0beta   <- ll_flexrsurv_beta0beta
    gr_G0A0B0AB    <- gr_ll_flexrsurv_GA0B0AB
    opgFunction <- opg_flexrsurv_G0A0B0AB
  }
  else {
    ll_G0A0B0AB    <- ll_flexrsurv_GA0B0AB_bh
    ll_gamma0      <- ll_flexrsurv_gamma0_bh
    ll_alpha0alpha <- ll_flexrsurv_alpha0alpha_bh
    ll_beta0beta   <- ll_flexrsurv_beta0beta_bh
    gr_G0A0B0AB    <- gr_ll_flexrsurv_GA0B0AB_bh
    opgFunction <- opg_flexrsurv_G0A0B0AB_bh
  }

  start <- list(gamma0= initgamma0, alpha0=initalpha0, beta0=initbeta0, alpha=initalpha, beta=initbeta)
  #initial value
  if(do.init){
    init_hessian <- FALSE
    # fit model to find gamma0, all aother parameters kept = 0
    if(do.init.gamma0){

      
      fit1 <- optim(par=initgamma0, fn=ll_gamma0, gr = NULL, 
                  method=method$optim_meth,
                  constraints=NULL, 
      # ll_flexrsurv_gamma0 args
                  alpha0=initalpha0, beta0=initbeta0, alpha=initalpha, beta=initbeta ,
                  Y=Y, X0=X0, X=X, Z=Z, 
                  expected_rate=expected_rate,
                      weights = weights,
                    step=step, Nstep=Nstep, 
                  intTD=intTD, intweightsfunc=intweightsfunc,
                  nT0basis=nT0basis,
                  Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                  nX0=nX0,
                  nX=nX, 
                  nTbasis=nTbasis,
                  Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
                  control = optim.control.gamma0, 
                  hessian = init_hessian, debug=debug.ll
                  )
      convgamma0 <- converged(fit1, step="during surch of init values for 'gamma0'")
      gamma0 <- fit1$par
      LL <- fit1$value
      # end of  if(do.init.gamma0)
    } else {
      gamma0 <- initgamma0
      LL <- +Inf
    }
    
      if(nalpha0+nalpha>0 & (do.init.alpha0 | do.init.alpha)) {
    # fit model to find alpha0 and alpha | previous gamm0 and beta0 and beta
        if (nZ) {
          initbeta <- rep(0,nbeta)
        }
        initalpha0alpha <- c(initalpha0, initalpha)
        fit2 <- optim(initalpha0alpha, fn=ll_alpha0alpha, gr = NULL, 
                    method=method$optim_meth,
                    constraints=NULL, 
    # ll_flexrsurv_alpha0alpha args
                    gamma0=gamma0, beta0=initbeta0, beta=initbeta,
                    Y=Y, X0=X0, X=X, Z=Z, 
                    expected_rate=expected_rate,
                      weights = weights,
                      step=step, Nstep=Nstep, 
                    intTD=intTD, intweightsfunc=intweightsfunc,
                    nT0basis=nT0basis,
                    Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                    ialpha0=Ialpha0, nX0=nX0,
                    nX=nX, 
                    ialpha= Ialpha,
                    nTbasis=nTbasis,
                    Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
                    control = optim.control.alpha,
                    debug=debug.ll,
                    hessian = init_hessian)

        convgalpha0alpha <- converged(fit2, step="during surch of init values for 'alpha'")
        LL <- fit2$value
      
        if(nX0){
          alpha0 <- fit2$par[Ialpha0]
        }
        if (nZ){
          alpha <- fit2$par[Ialpha]
        }
        # if initvalues for alpha found, update init value for beta
        if(do.init.alpha){
          do.init.beta <- TRUE
        }
      } else {
        alpha0 <- initalpha0
        alpha <- initalpha
        LL <- +Inf
      }
    
    
    if(nbeta0+nbeta>0 & (do.init.beta0 | do.init.beta)){
    # fit model with beta0 beta | gamma0 alpha0 alpha
        initbeta0beta <- c(initbeta0, initbeta)
        fit3<-optim(initbeta0beta, fn=ll_beta0beta, gr = NULL, 
                  method=method$optim_meth,
                  constraints=NULL, 
       # ll_flexrsurv_beta0beta args
                    alpha0=alpha0, alpha=alpha, gamma0=gamma0,
                    Y=Y, X0=X0, X=X, Z=Z, 
                    expected_rate=expected_rate,
                      weights = weights,
                    step=step, Nstep=Nstep, 
                    intTD=intTD, intweightsfunc=intweightsfunc,
                    nT0basis=nT0basis,
                    Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                    nX0=nX0,
                    ibeta0=Ibeta0, 
                    nX=nX, 
                    ibeta= Ibeta, 
                    nTbasis=nTbasis,
                    Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
                    control = optim.control.beta, 
                    debug=debug.ll,
                    hessian = init_hessian)
        convbeta0beta<-converged(fit3, step="during surch of init values for 'beta0' and 'bata'")
        if(nX){
          beta0<-fit3$par[Ibeta0]
        }
        if (nZ){
          beta <- fit3$par[Ibeta]
        }     
        LL <- fit3$value
      } else {
      beta0 <- initbeta0
      beta <- initbeta
      LL <- +Inf
    }
    
   # end   if(do.init)
  } else {
    gamma0 <- initgamma0
    alpha0 <- initalpha0 
    alpha  <- initalpha 
    beta0  <- initbeta0  
    beta   <- initbeta
    
    initG0A0B0AB <- c(gamma0, alpha0, beta0, alpha, beta)

    if(debug) cat("compute init LL value.\n")
    LL <- ll_G0A0B0AB(GA0B0AB= initG0A0B0AB, 
                      Y=Y, X0=X0, X=X, Z=Z, 
                      expected_rate=expected_rate,
                      weights = weights,
                      step=step, Nstep=Nstep, 
                      intTD=intTD, intweightsfunc=intweightsfunc,
                      nT0basis=nT0basis,
                      Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                      ialpha0=ialpha0, nX0=nX0,
                      ibeta0= ibeta0, nX=nX, 
                      ialpha=ialpha, 
                      ibeta= ibeta, 
                      nTbasis=nTbasis,
                      Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
                      debug=debug.ll
                      )
    if (debug) {
      cat("LL at init", LL, "\n")
	}
  }

   LLold<- LL
   conv <- TRUE
    iter <- -1
    algo <- "optim"
iter_hessian <- FALSE
optimfunction <- "optim"


    fit<-optim(initG0A0B0AB,
               fn=ll_G0A0B0AB,
               gr = gr_G0A0B0AB, 
               method=method$optim_meth,
               constraints=NULL, 
  # ll_GA0B0AB, args
               Y=Y, X0=X0, X=X, Z=Z, 
               expected_rate=expected_rate,
               weights=weights,
               step=step, Nstep=Nstep, 
               intTD=intTD, intweightsfunc=intweightsfunc,
               intTD_base=intTD_base,
               nT0basis=nT0basis,
               Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
               ialpha0=ialpha0, nX0=nX0,
               ibeta0= ibeta0, nX=nX, 
               ialpha=ialpha, 
               ibeta= ibeta, 
               nTbasis=nTbasis,
               Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
               control = optim.control, 
               debug=debug.ll,
               debug.gr=debug.gr,
                 hessian = dohessian)
  
  
    conv <- converged(fit, step="during maximisation")

    if(nT0basis){
      gamma0<-fit$par[1:nT0basis]
    }
    if (nX0){
      alpha0 <- fit$par[ialpha0]
    }
    if(nX){
      beta0<-fit$par[ibeta0]
    }
    if (nZ){
      alpha <- fit$par[ialpha]
      beta <- fit$par[ibeta]
    }
    LL <- fit$value
    iter <- NA
      
     
  # prepare returned value
    if(nT0basis){
      names(gamma0) <- names(fit$par[1:nT0basis])
    }

  if (nX0){
    names(alpha0) <- dimnames(X0)[[2]]
  }

  if(nX){
  namesbeta0 <- NULL
    for(ix in 1:nX) {
      namesbeta0 <- c(namesbeta0,
                      paste(paste("NPH(",
                                  dimnames(X)[[2]][ix],
                                  ", ",
                                  dimnames(Y)[[2]][1],
                                  "):",
                                  sep=""),
                            (2-Intercept_t_NPH[ix]):nTbasis,
                            sep="")
                      )
    }
    names(beta0) <-  namesbeta0
  }
    
  if (nZ){
    names(alpha) <-paste("NPHNLL:", dimnames(getDesignMatrix(Z))[[2]], sep="")
    namesbeta <- NULL
    for(iz in 1:nZ) {
      namesbeta <- c(namesbeta, 
                     paste(paste("NPHNLL(", getNames(Z), sep=""),
                            paste(dimnames(Y)[[2]][1],
                                  "):",
                                  2:nTbasis,
                                  sep=""),
                            sep="_"))
    names(beta) <-  namesbeta
    }
  }

   lcoef <- list(gamma0=gamma0, alpha0=alpha0, beta0=beta0 , alpha=alpha, beta=beta)
   coef <- c(lcoef$gamma0, lcoef$alpha0, lcoef$beta0 , lcoef$alpha, lcoef$beta)
 
# variance computation
    var <- NULL
    informationMatrix <- NULL
    cholinformationMatrix <- NULL
  if( vartype == "oim"){
# the variance matrix is obtained from the observed information matrix
# numericNHessian in package maxLik
    var <- NULL
    informationMatrix <- - fit$hessian
    options(show.error.messages = FALSE)
    cholinformationMatrix <- try(chol(informationMatrix), silent=TRUE)
    options(show.error.messages = TRUE)
    if( class(cholinformationMatrix)=="try-error"){
      var <- numeric(0)
      cat(geterrmessage())
    } else {
      options(show.error.messages = FALSE)
      var <- try( chol2inv(cholinformationMatrix) , silent=TRUE)
      options(show.error.messages = TRUE)
      if( class(var)=="try-error"){
        var <- numeric(0)
        cat(geterrmessage())
      }
            else{
        attr(var, "type") <- vartype
      }
    }
  } else if( vartype == "opg"){
    
    # the variance matrix is obtained from the outer product of the gradient
    # in the case of Type 1 censoring
    informationMatrix <- opgFunction(GA0B0AB=fit$par,
                                     Y=Y, X0=X0, X=X, Z=Z, 
                                     expected_rate=expected_rate,
                                     weights=weights,
                                     step=step, Nstep=Nstep, 
                                     intTD=intTD, intweightsfunc=intweightsfunc,
                                     intTD_base=intTD_base,
                                     nT0basis=nT0basis,
                                     Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                                     ialpha0=ialpha0, nX0=nX0,
                                     ibeta0= ibeta0, nX=nX, 
                                     ialpha=ialpha, 
                                     ibeta= ibeta, 
                                     nTbasis=nTbasis,
                                     Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
                                     debug.gr=debug.gr
                                     )
    var <- NULL
    cholinformationMatrix <- NULL
    options(show.error.messages = FALSE)
    cholinformationMatrix <- try(chol(informationMatrix), silent=TRUE)
    options(show.error.messages = TRUE)
    if( class(cholinformationMatrix)=="try-error"){
      cat(geterrmessage())
      var <- numeric(0)
    } else {
      options(show.error.messages = FALSE)
      var <- try( chol2inv(cholinformationMatrix) , silent=TRUE)
      options(show.error.messages = TRUE)
      if( class(var)=="try-error"){
        cat(geterrmessage())
        var <- numeric(0)
      }
      else {
        attr(var, "type") <- vartype
      }
    }
  } else {
    informationMatrix <- numeric(0)
    var <- numeric(0)
  }
   attr(informationMatrix, "type") <- vartype
   attr(var, "type") <- vartype

   
# if( !is.null(var)){
#   dimnames(var)[[1]] <-names(coef)
#   dimnames(var)[[2]] <-names(coef)
# }
  # computes the linear predictor objfit$linear.predictor
  linearPredictors <- .computeLinearPredictor_GA0B0AB(GA0B0AB=fit$par,
                                              Y=Y, X0=X0, X=X, Z=Z,
                                              nT0basis=nT0basis,
                                              Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                                              ialpha0=ialpha0, nX0=nX0,
                                              ibeta0= ibeta0, nX=nX, 
                                              ialpha=ialpha, 
                                              ibeta= ibeta, 
                                              nTbasis=nTbasis,
                                              Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
                                              debug=debug)
  
  # computes the fited rate
  fittedValues <- exp(linearPredictors)

# computes the cumulative rate??  
  cumulativeHazard <- .computeCumulativeHazard_GA0B0AB(GA0B0AB=fit$par,
                                             Y=Y, X0=X0, X=X, Z=Z,
                                             step=step, Nstep=Nstep, 
                                             intTD=intTD, intweightsfunc=intweightsfunc,
                                             nT0basis=nT0basis,
                                             Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                                             ialpha0=ialpha0, nX0=nX0,
                                             ibeta0= ibeta0, nX=nX, 
                                             ialpha=ialpha, 
                                             ibeta= ibeta, 
                                             nTbasis=nTbasis,
                                             Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
                                             debug=debug)

     # convergence assesment
     
  if (conv != TRUE) {
    warning("Ran out of iterations and did not converge")
  }
  
  res <- list(coefficients = coef,
              list_coefficients = lcoef,
              linear.predictors=linearPredictors,
              fitted.values=fittedValues,
              cumulative.hazard=cumulativeHazard, 
              var = var,
              informationMatrix=informationMatrix,
              loglik = LL,
              weights = weights,
              optimfunction=optimfunction,
              conv=conv,
              method = "flexrsurv.ll.fit",
              fit=fit,
              start=start)
  res
}

