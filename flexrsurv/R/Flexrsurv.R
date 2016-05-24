flexrsurv <- function(formula=formula(data),
                       data=parent.frame(), 
                       # knots.Bh and degree.Bh allow to define the Baseline Hazard
                       knots.Bh,
                       degree.Bh=3,
                       Spline=c("b-spline", "tp-spline", "tpi-spline"), # tp-spline for truncated power basis
                       log.Bh=FALSE,
                       bhlink=c("log", "identity"),
                       Min_T=0,
                       Max_T=NULL,
                       model=c("additive","multiplicative"),
                       rate=NULL, 
                       weights=NULL,
                       na.action=NULL,
                       int_meth=c("BANDS", "CAV_SIM", "SIM_3_8", "BOOLE"),
                       bands=NULL,
                       stept=NULL,              
                       init=NULL,
                       initbyglm=TRUE,
                       initbands=bands,
                       optim.control=list(trace=100, REPORT=1, fnscale=-1, maxit=25), 
                       optim_meth=c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "SANN", "Brent"),
                       control.glm=list(epsilon = 1e-8, maxit = 100, trace = FALSE, epsilon.glm = 1e-1, maxit.glm = 25 ),
                       vartype =  c("oim", "opg", "none"),
                       debug=FALSE
                       ){

  debug.ll <- FALSE
  debug.gr <- FALSE
  
  # formula: for example Surv(time,cens)~sex
  # data: the observed data set
  # rate: rate variable in data
    
  
  call   <- match.call() 
  ## extract x, y, etc from the model formula and frame
  if(missing(data)) data <- environment(formula)

  bhlink <- match.arg(bhlink) 
  vartype <- match.arg(vartype) 
  int_meth <- match.arg(int_meth) 
  optim_meth <- match.arg(optim_meth) 
  model  <- match.arg(model)     # type of model used (additive for remontet's model - multiplicative for mahboubi's model)
  Spline <- match.arg(Spline)    # choice of spline basis
                         
  # Validité des arguments
  #----------------------------------------------------------------------------------------

if (int_meth == "BANDS"){
	if( is.null(bands) ) {
		stop(gettextf("argument 'bands' must be specified if 'int_meth' = %s.", dQuote("BANDS"), domain=NA))
		}
	}
	else {
		if( is.null(stept) ) {
		stop(gettextf("argument 'step' must be specified if 'int_meth' != %s.", dQuote("BANDS"), domain=NA))
		}
	}
	
  
  mf <- match.call(expand.dots = FALSE)
  
  m <- match(c("formula", "data", "rate", "weights"), names(mf), 0L)
  if (m[1]==0)
    stop ("'formula' is required.")
  
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  specials <- c("NLL","NPH","NPHNLL", "nl", "td", "nltd")
  
  mf$formula <- if (missing(data)){ 
    terms(formula, specials=specials)
  }  else {
    terms(formula, data=data, specials=specials)
  }

  mf <- eval(mf, parent.frame())

  if (nrow(mf) == 0) 
    stop("No (non-missing) observations")    

  
  Y <- model.extract(mf, "response")
  if (!inherits(Y, "Surv")) 
    stop("Response must be a survival object.")
  
  Survtype <- attr(Y, "type")
  if ((ncol(Y) ==  2) && (Survtype != "right") ) {
    stop(gettextf("flexrsurv does not support %s type of censoring with (0, end] survival data", dQuote(Survtype), domaine=NA))
  } else if ((ncol(Y) ==  3) && (Survtype != "counting") ) {
      stop(gettextf("flexrsurv does not support %s type of censoring with (start, end] survival data", dQuote(Survtype), domaine=NA))
  }
  
  if (ncol(Y) ==  2)  {
    entervar <- rep(0, dim(Y)[1])   # enter time variable  
    timevar <- Y[, 1]   # time variable  
    fail    <- Y[, 2]   # dummy 1=event, 0=censur
  } else if (ncol(Y) ==  3){
    entervar <- Y[, 1]   # enter time variable  
    timevar <- Y[, 2]   # time variable  
    fail    <- Y[, 3]   # dummy 1=event, 0=censur
  }
  
  theweights <- as.vector(model.weights(mf))
  if (!is.null(theweights) && !is.numeric(theweights)) 
    stop("'weights' must be a numeric vector")
  if (!is.null(theweights) && any(theweights < 0)) 
    stop("negative weights not allowed")

  
  
  if (is.null(Min_T)) {
    Min_T <- max(entervar, min(bands) )
  }
  if (is.null(Max_T)) {
    Max_T <- max(c(bands, timevar))
  }

  therate <- as.vector(model.extract(mf, "rate"))
  if(!is.null(therate) && !is.numeric(therate))
      stop("'rate' must be a numeric vector")
  if (is.null(therate)){
    data$rate <- 0
  } else {
    data$rate <- therate
  }

  thedata <- data


  
  if (is.null(knots.Bh)) {
    stop("The knots for the baseline hazard must be specified in 'knots.Bh'.") 
  }

  #----------------------------------------------------------------------------------------
  # working formulas 
  #----------------------------------------------------------------------------------------

  # fix the formula
  # pass model and Spline parameter
  # check consistantcy of the formula
  # NLL(x) + NPH(X, t) or nl(x) +td(X,t)
  #----------------------------------------------------------------------------------------
  fixformula <- fix.flexrsurv.formula(formula=as.formula(formula), data=data, Spline=Spline, model=model, debug=debug)

  termsf <- terms(fixformula, specials=c("NLL","NPH","NPHNLL", "nl", "td", "nltd"))

  formulamle <- make.mle.formula(formula=fixformula, data= data)


      
  # Init values
  #----------------------------------------------------------------------------------------
  if (initbyglm == TRUE) {
    if (int_meth != "BANDS" & is.null(initbands)){
         initbands <- seq(Min_T, Max_T, stept*5)
		 warning("'initbands' has been set to seq(Min_T, Max_T, stept*5)", call. = FALSE)
      }

    
    name.runningtime <- ".t"
#  # build the formula for the GLM with the splited dataset
#  #----------------------------------------------------------------------------------------
    formulasplit <- make.glm.formula(formula=fixformula, data= data, name.runningtime=name.runningtime,
                                        Min_T=Min_T,Max_T=Max_T, model=model)
    splitdonnees <- split.data(jeudata=data, bands=initbands, fail=fail, entry=entervar, exit=timevar, name.runningtime=name.runningtime)
    if ((model=="additive") | ((model=="multiplicative") & (length(c(attr(termsf, "specials")$NPHNLL, attr(termsf, "specials")$nltp))==0)))  {
      control.glm <- do.call("control.flexrsurv.additive", control.glm)
      results <- flexrsurv.glm.fit(formula=formulasplit, data=splitdonnees, model=model,
                                Spline=Spline, degree.Bh=degree.Bh, knots.Bh=knots.Bh, log.Bh=log.Bh,
                                control=control.glm,
                                Min_T=Min_T, Max_T=Max_T, name.runningtime=name.runningtime, start=init)
           

    } else if((model=="multiplicative") &  (length(c(attr(termsf, "specials")$NPHNLL, attr(termsf, "specials")$nltp))!=0)) {

      control.glm <- do.call("control.flexrsurv.multiplicative", control.glm)
      results <- flexrsurv.glmiterative.fit(formula=formulasplit, data=splitdonnees, model=model, 
                                         Spline=Spline, degree.Bh=degree.Bh, knots.Bh=knots.Bh, log.Bh=log.Bh, 
                                         Min_T=Min_T, Max_T=Max_T, name.runningtime=name.runningtime, start=init,
                                         control=control.glm)

    }
    oldinit<-init
    init <- results$coef
    names(init) <- NULL
    if ( bhlink == "identity"){
      init[1:(length(knots.Bh)+degree.Bh+1)] <- exp(init[1:(length(knots.Bh)+degree.Bh+1)]) 
      }
    rm(results, splitdonnees)    
  } 
  # end of initbyglm

  # final Estimation,with variance computation
  #----------------------------------------------------------------------------------------
  call.ll <- match.call(definition=flexrsurv, expand.dots=FALSE)
  m <- match(names(formals(flexrsurv.ll)), names(call.ll), 0L)
  # keep flexrsurv.ll() args
  call.ll <- call.ll[c(1L, m)]
  call.ll[[1L]] <- quote(flexrsurv::flexrsurv.ll)
  call.ll$formula <- formulamle
  if (initbyglm == TRUE) {
  # rem: there is stil a pb in case of factor variables in the linear effects
    call.ll$init <- init[!is.na(init)]
 #   call.ll$init <- init
  }
  results <- eval(call.ll, parent.frame()) 


    if( debug == FALSE ){
      results$fit <- NULL
    }
    
    
  
      

  # Affichage de la formule initiale et de la formule
  #----------------------------------------------------------------------------------------
  results[["call"]] <- call
  results[["formula"]] <- fixformula
  if(dim(Y)[2]==3){
    results[["entertime"]] <- entervar
  }
  results[["time"]] <- timevar
  results[["workingformula"]] <- formulamle
  

  
  # Simplified names of coefficients
  #----------------------------------------------------------------------------------------

  names(results$coefficients) <- make.shortnames.coefficients(names(results$coefficients),
                                                              formula=results[["workingformula"]] ,
                                                              model=model,
                                                              Spline=Spline,
                                                              knots.Bh=knots.Bh,
                                                              degree.Bh=degree.Bh,
                                                              log.Bh=log.Bh)


  
  if( !is.null(dim(results$var))){
    dimnames(results$var)[[1]] <-  names(results$coefficients)
    dimnames(results$var)[[2]] <-  names(results$coefficients)
  }
  if( !is.null(dim(results$informationMatrix))){
    dimnames(results$informationMatrix)[[1]] <-  names(results$coefficients)
    dimnames(results$informationMatrix)[[2]] <-  names(results$coefficients)
  }
  
  results$data <- thedata
  results$rate <- therate

  results$terms <- terms(formulamle, special=c("NLL", "NPH", "NPHNLL", "td", "nltd"),
                         data = data)

  if (initbyglm==TRUE) {
      results$init<- oldinit
      results$control.glm <- control.glm
  } else {
      results$init<- init
  }
  
  class(results) <- c("flexrsurv", class(results))
  
  return(results)
  
}



