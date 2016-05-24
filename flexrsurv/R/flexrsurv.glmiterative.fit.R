# version du 15 oct 2014

# Function of estimation using the iterative alternative algorithm, in case of NPHNLL effects
#----------------------------------------------------------------------------------------

flexrsurv.glmiterative.fit <- function(formula, data,
                                    Spline,
                                    degree.Bh=3,
                                    knots.Bh,
                                    log.Bh=FALSE, 
                                    model,
                                    control=list(epsilon = 1e-8, maxit = 100, trace = FALSE, epsilon.glm = 1e-2, maxit.glm=NULL),
                                    Min_T, Max_T, 
                                    name.runningtime=".t",
                                    start=NULL) {

  
  tik <- data$tik
  data$exp_nbevent <-  tik*data$rate

  special <- c("NPH","NLL", "NPHNLL") 
  Terms <- if (missing(data)){
    terms(formula, special, keep.order = TRUE)
  } else {
    terms(formula, special, data = data, keep.order = TRUE)
  }

  
  theTerms <- attr(Terms, "term.labels")
  n <- length(theTerms)                               # number of effect (i=1:n)

  
  # tolerance for the alternative algorithm
  epsilon.alternative <- control$epsilon
  maxit.alternative <- control$maxit

  trace <- control$trace
  
  # tolerance for the inner glm
  control.glm <- control
  control.glm$epsilon <- control$epsilon.glm
  control.glm$maxit   <- max(control$maxit.glm ,1)
  control.glm$epsilon.glm <- NULL
  control.glm$maxit.glm   <- NULL
  control.glm <- do.call("glm.control", control.glm)
  
  # Definition of the link function
  #----------------------------------------------------------------------------------------
  # link function
  linkfun  <- function(mu){
    log(mu - data$exp_nbevent)
  }
  # inverse of the link function
  linkinv  <- function(eta){
    pmax(exp(eta), .Machine$double.eps/2) + data$exp_nbevent 
  }
  # dmu/deta
  mu.eta   <- function(eta){
    pmax(exp(eta), .Machine$double.eps/2)
  }
  valideta <- function(eta){
    TRUE
  }
  flexrsurv.link <- structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
                              valideta = valideta, name = "flexrsurv"), class = "link-glm")
  
  # Estimation
  # At least one NPHNLL effect tested --> iterative alternative algorithm (alpha, beta, ...)
  #--------------------------------------------------------------------------------------
  
  # if none NPHNLL effect, use estimation.glm()

  index.NPHNLL <- grep("NPHNLL", theTerms) 
  is.NPHNLL <- rep(FALSE, n)
  is.NPHNLL[index.NPHNLL] <- TRUE
  
  # NPHNLL(x,t) =alpha(x)beta(t)
  # alpha(x) = NLL(x)
  # beta(t)  = b1(t) + NLL(T, intercept=FALSE)
  #      where b1 is the veryfirst spline basis
  #            if tp-spline b1(t) = 1 ,
  #            if b-spline b1(t) is the only basis such that b1(left_boundary_knot) != 0   
  # 
  # NPH step : alpha are assumed known, estimation of beta
  # in NPH step, NPHNLL(x, t) is NPH(Alpha(x), t)
  # NLL step   : beta are assumed known, estimate alpha
  # in NLL step, NPHNLL(x, t) is  NLLbeta(beta(t), x)=NLL(x)*beta(t)
  # initialization  

  names_alphax <- vector("character", n)      # alpha(x) variable 
  names_alphaxb1 <- vector("character", n)    # alpha(x) variable *1st spline basis
  names_betaT <- vector("character", n)       # beta(t) variable 
  
  nbNPH <- vector("double",n)                       # nbNPH[i] : nb of coef of the ith effect in step NPH 
  nbNLL <- vector("double",n)                       # nbNLL[i] : nb of coef of the ith effect in step NLL 
  nbcoeffeffect <- vector("double",n)             # nbcoeffeffect [i] : total nb of coeff of the ith effect in the model

  SnbNPH <- vector("double",n)                      # SnbNPH[i] : 1- index of the 1st coef of the ith effect estimated at step NPH
  SnbNLL <- vector("double",n)                      # SnbNLL[i] : 1- index of the 1st coef of the ith effect estimated at step NLL
  Snbcoeffeffect <- vector("double",n)            # Snbcoeffeffect [i] : 1- index of the 1st coef of the ith effect estimated in the model
  
  # X[[i]] & TT[[i]]: design matrix corresponding to the ith effect
  XX <- vector("list",n) ;
  TT <- vector("list",n) ;
  
  # alpha[[i]] = valeurs d'initialisation des coefficients alpha pour le ième effet, ...
  alpha <- vector("list",n)                       # if is.NPHNLL[i], alpha[[i]] is the running set of coef of alpha(x) 
  beta <- vector("list",n)                        # if is.NPHNLL[i], beta[[i]]  is the running set of coef of beta(t)
  ALPHAx <- matrix(0,ncol=n, nrow=dim(data)[1])   # if is.NPHNLL[i], ALPHAx[[i]] = alpha(x) = XX %*% alpha + min(x)
  BETAt <- matrix(0,ncol=n, nrow=dim(data)[1])    # if is.NPHNLL[i], BETAt[[i]]  = beta(t)  = TT %*% c(1, beta)
  
  #----------------------------------------------------------------------------------------
  
  # Definition of the baseline : gamma0_t
  #----------------------------------------------------------------------------------------

  if (length(knots.Bh)>=2) {
    k<-(knots.Bh)
    for (l in 2:length(k)) {
      k[l] <- paste(k[l-1],k[l],sep=",") 
    } 
    knots.Bh <- paste("c(",k[length(k)],")",sep="")  # consider multiple knots
  }
  
  gamma0_t <- as.character(paste("NLL(", name.runningtime , ', Spline="', Spline,'", Knots=', knots.Bh,
                                 ", Degree=", degree.Bh, ", Log=", log.Bh,
                                 ", Intercept=TRUE, Boundary.knots=c(", Min_T, ", ", Max_T, "))", sep=""))
  
  
  # get index/number of coef of the baseline hazard
  #---------------------------------------------------------------------------------------- 
  
  f_Bh    <- paste(".fail ~ ", gamma0_t, "-1", sep = "") 
  SnbNLL[1] <- dim(design(f_Bh,data))[2]
  SnbNPH[1] <- SnbNLL[1]
  Snbcoeffeffect[1] <- SnbNLL[1]


  #----------------------------------------------------------------------------------------
  #  formulas used in the alternating conditional algorithm
  #----------------------------------------------------------------------------------------

  formulNPH <- make.formulastepNPH.terms(Terms, data, baseline=gamma0_t, response=".fail", tik="tik",)
  formulNLL <- make.formulastepNLL.terms(Terms, data, baseline=gamma0_t, response=".fail", tik="tik",)

  # get pline parameters of NPHNLL effects
  Degree <- vector("double",n)
  Degree.t <- vector("double",n)
  Boundary.Knots <- vector("list",n)   
  Knots <- vector("list",n)   
  Intercept <- vector("logical",n)   
  Knots.t <- vector("list",n)
  Boundary.Knots.t <- vector("list",n)   
  Intercept.t <- vector("logical",n)   
  xnames <- vector("character",n)                 # xnames[i] name of the variable of the ith efect in the formula
  tnames <- vector("character",n)                 # tnames[i] name of the time variable of the ith efect in the formula
  coefnames <- vector("list", n)

  if(length(index.NPHNLL) >0){
    for (i in attr(Terms, "specials")["NPHNLL"]){
        thecall <-  match.call(NPHNLL, attr(Terms,"variables")[[i+1]])
        indxterm <- variable2term(i, Terms)
        xname <- as.character(thecall[[2]])
        tname <- as.character(thecall[[3]])
        tnames[i] <- tname
        valDegree <- eval(as.expression(thecall[["Degree"]]))
        if (is.null(valDegree)) {
          valDegree <- 3
        }
        Degree[indxterm] <- valDegree            
        valDegree.t <- eval(as.expression(thecall[["Degree.t"]]))
        if (is.null(valDegree.t)) {
          valDegree.t <- 3
        }
        Degree.t[indxterm] <- valDegree.t          
        if (!is.null(eval(as.expression(thecall[["Knots"]])))) {
          Knots[[indxterm]] <- eval(as.expression(thecall[["Knots"]]))
        }
        if (!is.null(eval(as.expression(thecall[["Intercept"]])))) {
          Intercept[indxterm] <- eval(as.expression(thecall[["Intercept"]]))
        }
        else {
         Intercept[indxterm] <- FALSE
        }
        if (!is.null(eval(as.expression(thecall[["Boundary.knots"]])))) {
          Boundary.Knots[[indxterm]] <- eval(as.expression(thecall[["Boundary.knots"]]))
        }
        else {
           Boundary.Knots[[indxterm]] <-  with(data, eval(call("range", as.name(xname))))
        }
        
        if (!is.null(eval(as.expression(thecall[["Knots.t"]])))) {
          Knots.t[[indxterm]] <- eval(as.expression(thecall[["Knots.t"]]))
        }
        if (!is.null(eval(as.expression(thecall[["Intercept.t"]])))) {
          Intercept.t[indxterm] <- eval(as.expression(thecall[["Intercept.t"]]))
        }
        else {
          Intercept.t[indxterm] <- FALSE
        }
        Boundary.Knots.t[[indxterm]] <- range( Min_T, Max_T)
        if (!is.null(eval(as.expression(thecall[["Boundary.knots.t"]])))) {
          Boundary.Knots.t[[indxterm]] <- eval(as.expression(thecall[["Boundary.knots.t"]]))
        }
        coefnames[[indxterm]] <- c(paste(paste("NPHNLL(", xname, ", ", tname, ")", xname, ":", sep=""),
                                         1:(length(Knots[[indxterm]]) + Degree[indxterm] + Intercept[indxterm]), sep=""),
                                   paste(paste("NPHNLL(", xname, ", ", tname, ")", tname, ":", sep=""),
                                         1:(length(Knots.t[[indxterm]]) + Degree.t[indxterm] + Intercept.t[indxterm]), sep=""))
        rm(thecall,valDegree,valDegree.t, xname, tname)                          
      }
  }

  #----------------------------------------------------------------------------------------
  # get the NAMES OF the VARIABLES 
  #----------------------------------------------------------------------------------------
  for (i in (3:length(attr(Terms, "variables")))){
    indxterm <- variable2term(i-1, Terms)
    if(is.name(attr(Terms, "variables")[[i]])){
      for( j in indxterm){
        if(xnames[j] == ""){
          xnames[j] <- as.character(attr(Terms, "variables")[[i]])
        }
        else {
          xnames[j] <- paste(xnames[j], ":", as.character(attr(Terms, "variables")[[i]]), sep="")
        }
      }
    }
    else if(is.call(attr(Terms, "variables")[[i]])){
      fun <- match.fun(attr(Terms, "variables")[[i]][[1]])
      thecall <-  match.call(fun, attr(Terms,"variables")[[i]])
      for( j in indxterm){
        if(xnames[j] == ""){
          xnames[j] <- as.character(thecall)[[2]]
        }
        else {
          xnames[j] <- paste(xnames[j], ":", as.character(thecall)[[2]], sep="")
        }
      }
    }
    else {
      stop("the formula is not correct")
    }
  }

  for (i in 1:n) {
    if (is.NPHNLL[i]) {
      names_alphax[i] <- paste("alpha",xnames[i],sep="")                                    
      names_alphaxb1[i] <- paste("alpha",xnames[i],"b1",sep="")                                    
      names_betaT[i] <- paste("betaT",xnames[i],sep="")
    }
  }


  #--------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  #-------- build index/number of coef of each effect for the 2 steps and for the total model -------------------------------
  #---------and  get the spline basis (X and T) of the NPHNLL effects
  #--------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  # vector of names of variables, nbNPH, nbNLL, SnbNPH, SnbNLL (use to define start.NPH in estim_stepNPH)
  #----------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------

  formulNPH2 <- NPHNLL2NPHalpha.formula(formula, data, response=".fail")
  formulNLL2 <- NPHNLL2NLL.formula(formula, data, response=".fail")

  TermsNPH <- if (missing(data)){
    terms(formulNPH2, special, keep.order = TRUE)
  } else {
    terms(formulNPH2, special, data = data, keep.order = TRUE)
  }
  theTermsNPH <- attr(TermsNPH, "term.labels")
  
  TermsNLL <- if (missing(data)){
    terms(formulNLL2, special, keep.order = TRUE)
  } else {
    terms(formulNLL2, special, data = data, keep.order = TRUE)
  }
  theTermsNLL <- attr(TermsNLL, "term.labels")

  
  nb_coef_beta <- 0L
  for (i in 1:n) {
    tmpformula <- paste(".fail ~ ", theTermsNLL[i], "-1", sep = "")
    designmat <- design(tmpformula,data)
    nbNLL[i] <- dim(designmat)[[2]]
    if (is.NPHNLL[i]) {
      # get XX design matrix of NPHNLL effects to compute Alpha(X)
      XX[[i]] <- designmat
#      dimnames(XX[[i]]) <- NULL

      # get TT design matrix of NPHNLL effects to compute beta(t)
      # TT is got from the matrix design if formulaNPH + first basis (ie intercept.t=TRUE)
      tmpformula <- paste(".fail ~ ", theTermsNPH[i], "-1", sep = "")
      TT[[i]] <- design(tmpformula,data)  

#      dimnames(TT[[i]]) <- NULL
      # remove intercept (the first asis is absent in the main model
      nbNPH[i] <- dim(TT[[i]])[[2]]-1
      nb_coef_beta <- nb_coef_beta + dim(TT[[i]])[[2]]-1
      nbcoeffeffect[i] <- nbNPH[i] + nbNLL[i]
    }
    else {
        # else same as nbNLL
      nbNPH[i] <- nbNLL[i]
      nbcoeffeffect[i] <- nbNPH[i] 
    }
        # Snb???[i] : index of the 1st coeff (estimated) corresponding to the ith effect
    if(i > 1){
      SnbNLL[i]   <- SnbNLL[i-1]+nbNLL[i-1]
      SnbNPH[i] <- SnbNPH[i-1]+nbNPH[i-1]
      Snbcoeffeffect[i] <- Snbcoeffeffect[i-1]+nbcoeffeffect[i-1]
    }
  }  # for i in 1:n
  ncoeff <- sum(nbcoeffeffect) + SnbNLL[1]        
  ncoeff.stepNPH <- sum(nbNPH) + SnbNLL[1]        
  ncoeff.stepNLL  <- sum(nbNLL) + SnbNLL[1]        
  
  
  #--------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  #-------- build index/number of coef of each effect for the 2 steps and for the total model -------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------------------------
  # vector of names of variables, nbNPH, nbNLL, SnbNPH, SnbNLL (use to define start.NPH in estim_stepNPH)
    
  
      #----------------------------------------------------------------------------------------
      # starting values of coefficients
      #----------------------------------------------------------------------------------------
  if (!is.null(start)) {
    # dispatch starting values in start.stepNPH end start.stepNLL
    # and in coef.stepNPH end coef.stepNLL
        # control of length of vector "start"  
    
    start.NPH <- rep(0, ncoeff.stepNPH)
    start.NLL <- rep(0, ncoeff.stepNLL)
    
    if (length(start)!=ncoeff) {
      stop(gettextf("Wrong length for inital values. Length of init values must be equal to %i.\n", ncoeff, domaine=NA))
    }
        
        # define start.NPH
      #baseline gamma0(t)
    start.NPH[1:(SnbNPH[1])] <- start[1:(SnbNPH[1])]
    start.NLL[1:(SnbNLL[1])]  <- start[1:(SnbNLL[1])]

    for (i in 1:n) {
      if (is.NPHNLL[i]) {
              #in coefficients, fires alpha(x) then betat(t)
              # get beta
        beta[[i]] <- start[Snbcoeffeffect[i]+nbNLL[i]+(1:nbNPH[i])]
        start.NPH[SnbNPH[i]+(1:nbNPH[i])] <- start[Snbcoeffeffect[i]+nbNLL[i]+(1:nbNPH[i])]
              # get alpha
        alpha[[i]] <- start[Snbcoeffeffect[i]+(1:nbNLL[i])]
        start.NLL[SnbNLL[i]+(1:nbNLL[i])] <- start[Snbcoeffeffect[i]+(1:nbNLL[i])]
      }
      else { # !NPHNLL
        start.NPH[SnbNPH[i]+(1:nbNPH[i])] <- start[Snbcoeffeffect[i]+(1:nbNPH[i])]
        start.NLL[SnbNLL[i]+(1:nbNLL[i])]  <- start[Snbcoeffeffect[i]+(1:nbNLL[i])]
      } 
    } 
  } # end   if (!is.null(start)) 
  else {
    # no starting values supplied
    start.NPH <- NULL
    start.NLL <- rep(0, ncoeff.stepNLL)

    # But alpha(x) must be initialised 
    for (i in 1:n) {
      if (is.NPHNLL[i]) {
      #  NPHNLL(x, t) is assumed to be X (ie alpha(x)=x, beta(t) =1 + 0(t))
        if (Spline=="tp-spline") {
          alpha[[i]] <- c(1,rep(0, nbNLL[i]-1))
          beta[[i]] <- rep(0, nbNPH[i])
        }
        else {
          # when no starting values, alpha(X) is assumed equal to X-Xmin, that is Alpha(Xmin)=0 and alpha(x) linear),
          # where Xmin is the lowest boundary.knot
          # thus using Marsden's identities, alpha[[i]] = marsdens.matrix[[i]] %*% Knots[[i]]
              
          Knots.marsden <- sort(c(rep(Boundary.Knots[[i]], Degree[i]+1), Knots[[i]]))
          marsden.matrix <- matrix(0, ncol=length(Knots[[i]]) + 2*(Degree[i]+1), nrow=length(Knots[[i]])+Degree[i]+1)
          
          for (j in 1:(length(Knots[[i]])+Degree[i]+1)) {
            marsden.matrix[j, j+(1:Degree[i])] <- 1.0
          }

          tmpval <- marsden.matrix%*% Knots.marsden/Degree[i] - Boundary.Knots[[i]][1]
          alpha[[i]] <- tmpval[-1]
        } # start = null
        start.NLL[SnbNLL[i]+(1:nbNLL[i])] <- alpha[[i]]
      } # effet NPHNLL 
    } # end of for (i in 1:n) # 
  } # else
      
      

  
  # Iterative alternative algorithm : 
  # (2 steps :
  #     -1-(step NPH)- estimation of beta (=with alpha known), 
  #     -2-(step NLL) then estimation of alpha (=with beta known)) 
  #----------------------------------------------------------------------------------------------------------------------------
  # Initialization of the log-likelihood
  #----------------------------------------------------------------------------------------
  
  ll    <- +Inf
  llold <- +Inf
  conv <- FALSE

  for( iter in 1L:maxit.alternative) {   

    
    #--------------------#
    # STEP NPH: ALPHA KNOWN #
    #--------------------#
    
    
    
    # update ALPHA(X) using previous estimate
    #--------------------------------------------------------
    
    for (i in 1:n) {
      if (is.NPHNLL[i]) { 
        if (Spline=="b-spline") {
#          ALPHAx[,i] <- XX[[i]]%*%alpha[[i]] + Boundary.Knots[[i]][1]
          ALPHAx[,i] <- XX[[i]]%*%alpha[[i]] 
        }
        else {
          ALPHAx[,i] <- XX[[i]]%*%alpha[[i]]
        }
        data[,names_alphax[i]] <- ALPHAx[,i]           

        # offset for NPHNLL(x, T) is alpha(x)*b1(T)
        if (Spline=="b-spline") {
          data[,names_alphaxb1[i]] <- ALPHAx[,i]*TT[[i]][,1] 
        }
        else {
          # because b1(T) =1
          data[,names_alphaxb1[i]] <- ALPHAx[,i]
       } 
#        cat("names data", names(data), "\n", sep=" ")
#        save(data, file=paste("iterationglm", runif(1), ".RData", sep=""))
          
      } # NPHNLL 
    } # end of loop for (number of effect), computation of ALPHA(X)
    
    
    
    # update starting values of the glm for NPH step
    #--------------------------------------------------------    
    
    if (iter >1) { # start at other iteration
      start.NPH <- estim_stepNPH$coeff
      
      #update baseline gamma0(t)
      start.NPH[1:(SnbNPH[1])] <- estim_stepNLL$coeff[1:(SnbNLL[1])]
        
       #keep coefficients of NPHNLL, update the other
      for (i in 1:n) {
        if (!is.NPHNLL[i]) {
          start.NPH[SnbNPH[i]+(1:nbNPH[i])] <- estim_stepNLL$coeff[SnbNLL[i]+(1:nbNLL[i])]
        }
      }
    } # iter >1


    estim_stepNPH <- glm(formula=formulNPH, family=poisson(link=flexrsurv.link), data=data,
                          control=control.glm, start=start.NPH)
    expetahat <- pmax(exp(predict(estim_stepNPH)), .Machine$double.eps)
    
    llinterm <- sum(-expetahat + data$.fail*log(expetahat/tik + data$rate))

    if (trace == TRUE ){
      cat("iteration ", iter, ": ll=", llinterm, "\n")
    }
      # Save the ubdated coefficients beta  (under contraint beta0=1)
    #----------------------------------------------------------------------------------------
    #- put aliased coeff to 0
    for (j in 1:length(estim_stepNPH$coeff)) {
      if (is.na(estim_stepNPH$coeff[j])) {
        estim_stepNPH$coeff[j]<-0
      }
    }

    for (i in 1:n) {                    
      if (is.NPHNLL[i]) { 
          beta[[i]] <- c(1, estim_stepNPH$coeff[SnbNPH[i]+(1:nbNPH[i])])
        }
      else {
        beta[[i]] <- estim_stepNPH$coeff[SnbNPH[i]+(1:nbNPH[i])]
      }        
    }
    #-------------------#
    # STEP2: NLL step,  BETA KNOWN #
    #-------------------#   

  # update beta(T)
  for (i in 1:n) {
    if (is.NPHNLL[i]) { 
      # BETAtx : new variable to be substituted in x (when beta are known)
      BETAt[,i] <- (TT[[i]]%*%beta[[i]])
      data[,names_betaT[i]] <- BETAt[,i]
    }
  } 
    
    # update starting values of the glm for NLL step
    #--------------------------------------------------------    
    
    if (iter >1) { 
      start.NLL <- estim_stepNLL$coeff
    }
      #update baseline gamma0(t)
    start.NLL[1:(SnbNLL[1])] <- estim_stepNPH$coeff[1:(SnbNPH[1])]
        
       #keep coefficients of NPHNLL, update the other
    for (i in 1:n) {
      if (!is.NPHNLL[i]) {
        start.NLL[SnbNLL[i]+(1:nbNLL[i])] <- estim_stepNPH$coeff[SnbNPH[i]+(1:nbNPH[i])]
      }
    }
      

    
    estim_stepNLL <- glm(formula=formulNLL, family=poisson(link=flexrsurv.link), data=data, 
                         control=control.glm, start=start.NLL)
    
    
    # Save the coefficients alpha
    #----------------------------------------------------------------------------------------
    #- put aliased coeff to 0
    for (j in 1:length(estim_stepNLL$coeff)) {
      if (is.na(estim_stepNLL$coeff[j])) {
        estim_stepNLL$coeff[j]<-0
      }
    }  
    
    for (i in 1:n) {                    
      if (is.NPHNLL[i]) { 
          alpha[[i]] <- estim_stepNLL$coeff[SnbNLL[i]+(1:nbNLL[i])]
        }
      else {
        alpha[[i]] <- estim_stepNLL$coeff[SnbNLL[i]+(1:nbNLL[i])]
      }        
    }

    # compute log likelyhood
    expetahat <- pmax(exp(predict(estim_stepNLL)), .Machine$double.eps)

    ll <- sum(-expetahat + data$.fail*log(expetahat/tik + data$rate))

    if (trace == TRUE ){
      cat("iteration ", iter, ": ll=", ll, "\n")
    }

    if (abs(ll - llold)/(0.1 + abs(ll)) < epsilon.alternative ) {
      conv <- TRUE
      break
    }
    else {
      llold <- ll
    }

  } # for(iter ...)

  #-------------------------------#
  # convergence assessment        #  
  #-------------------------------#
  

  if (!conv){
    warning("algorithm did not converge", call. = TRUE)
  }
    
  
  #-------------------------------#
  # Fusion of the 2 glm             #  
  #-------------------------------#
  
  # get information from the last NLL step
  estim <- estim_stepNLL

  # buid the whole set of coeff
  estim$coefficients <- rep(0, ncoeff)
  #baseline gamma(0)
  estim$coefficients[1:SnbNLL[1]] <- estim_stepNLL$coefficients[1:SnbNLL[1]]
  names(estim$coefficients)[1:SnbNLL[1]] <- names(estim_stepNLL$coefficients)[1:SnbNLL[1]]
  for (i in 1:n) {
    if (is.NPHNLL[i]) {
              #in coefficients, first alpha(x) then betat(t)
              # get alpha
      estim$coefficients[Snbcoeffeffect[i]+(1:nbNLL[i])] <- estim_stepNLL$coefficients[SnbNLL[i]+(1:nbNLL[i])] 
              # get beta
      estim$coefficients[Snbcoeffeffect[i]+nbNLL[i]+(1:nbNPH[i])] <- estim_stepNPH$coefficients[SnbNPH[i]+(1:nbNPH[i])] 
      names(estim$coefficients)[Snbcoeffeffect[i]+(1:(nbNPH[i]+nbNLL[i]))] <- coefnames[[i]]
    }
    else { # !NPHNLL
      estim$coefficients[Snbcoeffeffect[i]+(1:nbNLL[i])] <- estim_stepNLL$coefficients[SnbNLL[i]+(1:nbNLL[i])] 
      names(estim$coefficients)[Snbcoeffeffect[i]+(1:nbNLL[i])] <- names(estim_stepNLL$coefficients)[SnbNLL[i]+(1:nbNLL[i])] 
    } 
  }
  
  # Others elements 
  #----------------------------------------------------------------------------------------          
  
  # number of made iteration
  estim$iter <- iter
  estim$formula <- formula  
  
  # offset considered
  estim$offset <- (estim_stepNPH$offset+estim_stepNLL$offset)-log(data$tik) 
  estim$loglik <- ll
  estim$workingformulaNLL <- formulNLL
  estim$workingformulaNPH <- formulNPH

  estim$rank <- estim$rank + nb_coef_beta

  
  estim$residuals <- NULL
  estim$fitted.values <- NULL
  estim$linear.predictors <- NULL
  estim$deviance <- NULL
  estim$aic <- NULL
  estim$null.deviance <- NULL
  estim$weights <- NULL
  estim$prior.weights <- NULL
  estim$df.null <- NULL
  estim$offset <- NULL


  estim$effects <- NULL 
  estim$R <- NULL
  estim$qr <- NULL
  estim$model <- NULL
  estim$coeff <- NULL
  


  
  class(estim) <- c("glm.iterative")
  
  
  return(estim)     
} 
