flexrsurv.ll <- function(formula=formula(data), 
                         data=parent.frame(), 
                         # knots.Bh and degree.Bh allow to define the Baseline Hazard
                         knots.Bh=NULL,   
                         degree.Bh=3,
                         Spline=c("b-spline","tp-spline","tpi-spline"), # tp-spline for truncated power basis
                         log.Bh=FALSE,
                         bhlink=c("log", "identity"),
                         Min_T=0,
                         Max_T=NULL,
                         model=c("additive","multiplicative"),
                         rate=NULL, 
                         weights=NULL,
                         na.action=NULL, 
                         int_meth=c("CAV_SIM", "SIM_3_8", "BOOLE", "GLM", "BANDS"),
                         stept=NULL,              
                         bands=NULL,
                         init=NULL,
                         optim.control=list(trace=100, REPORT=1, fnscale=-1, maxit=25), 
                         optim_meth=c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "SANN", "Brent"),
                         vartype =  c("oim", "opg", "none"),
                         debug=FALSE
                         ){

debug.gr <- debug.ll <- FALSE
  if( debug > 100000){
  debug.gr <- debug%%1000 + debug%/%100000 * 1000 
}
if( debug > 10000){
  debug.ll <- debug%%1000 + (debug%/%10000)%%10 * 1000 
}
  debug <- debug%%1000 + (debug%/%1000)%%10 * 1000

# force optim.control$fnscale to -1 (maximum of the objective function in optim)
  optim.control$fnscale <- -1
  vartype  <- match.arg(vartype)       # type of variance
  optim_meth <- match.arg(optim_meth)
  int_meth <- match.arg(int_meth)
  if (int_meth == "BANDS" ){
	int_meth = "GLM"
  }
# setting up variables
    call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  
  m <- match(c("formula", "data", "rate", "weights"), names(mf), 0L)
  if (m[1]==0)
    stop ("The formula argument is required.")
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
    
  special <- c("NPH","NLL", "NPHNLL", "nl", "td", "nltd") 
  Terms <- if (missing(data)){
    terms(formula, specials=special)
  } else {
    terms(formula, specials=special, data = data)
  }
  
  mf$formula <- Terms


  mf <- eval(mf,  sys.parent())
  mt <- attr(mf, "terms")
 
  # remove data containing NAs
  na.act <- attr(m, "na.action")
  if( !is.null( na.act )) {
    data <- data[na.act,]
  }
  
  
    # pour l'intercept, voir avec l'option utilisée dans gamma0(t)
    intercept <- attr(mt, "intercept")

    Y <- model.extract(mf, "response")
    if (!inherits(Y, "Surv")) {
      stop("Response must be a survival object")
    }
    type <- attr(Y, "type")
  
    if ((ncol(Y) ==  2) && (type != "right") ) {
      stop(gettextf("flexrsurv does not support %s type of censoring with (0, end] survival data", dQuote(type), domaine=NA))
    } else if ((ncol(Y) ==  3) && (type != "counting") ) {
      stop(gettextf("flexrsurv does not support %s type of censoring with (start, end] survival data", dQuote(type), domaine=NA))
    }
        
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights)) 
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0)) 
    stop("negative weights not allowed")

  rate <- as.vector(model.extract(mf, "rate"))
  if(!is.null(rate) && !is.numeric(rate))
      stop("'rate' must be a numeric vector")
  if (is.null(rate)){
    data$rate <- 0
  } else {
    data$rate <- rate
  }



# build Design Matrics

# get T-splines

      if(is.null(Max_T)){
        Max_T <- max(Y[,1:(ncol(Y)-1)])
      }
      if(is.null(Min_T)){
        Min_T <- 0
      }

if(Spline=="b-spline"){
  
  Spline_t0 <- MSplineBasis(knots=c(Min_T, knots.Bh, Max_T),
                            degree=degree.Bh,
                            keep.duplicates=TRUE,
                            log=log.Bh)

  # no log basis for NPH and NPHNLL and td and nltd effects
  Spline_t<- MSplineBasis(knots=c(Min_T, knots.Bh, Max_T),
                          degree=degree.Bh,
                          keep.duplicates=TRUE,
                            log=FALSE)
} else if(Spline=="tp-spline" || Spline=="tpi-spline"){
# because time is >0, tp-splines are naturaly increasing
      Spline_t0 <-TPSplineBasis(knots=knots.Bh,
                                degree=degree.Bh,
                                min=Min_T,
                                max=Max_T,
                                log=log.Bh,
                                type="standard")

  # no log basis for NPH and NPHNLL  and td and nltd effects
      Spline_t <-TPSplineBasis(knots=knots.Bh,
                               degree=degree.Bh,
                                min=Min_T,
                                max=Max_T,
                               log=FALSE,
                                type="standard")
    } 

# lists of variables in the model/formula
    list_var_LIN <- all_LIN_vars(Terms)
    list_var_NLL <- all_specials_vars(Terms, specials="NLL",
                           unique = FALSE,
                           order="formula")
    list_var_NPH <- all_specials_vars(Terms, specials="NPH",
                           unique = FALSE,
                           order="formula")
    list_var_NPHNLL <- all_specials_vars(Terms, specials="NPHNLL",
                           unique = FALSE,
                           order="formula")
    
# coherence between LIN, NLL, NPH, NPHNLL terms
      
    var_LIN_NPHNLL <- list_var_LIN %in% c(list_var_NLL, list_var_NPH, list_var_NPHNLL)
    var_NPH_NLL <- list_var_NPH %in% list_var_NLL
    var_NLL_NPHNLL <- list_var_NLL %in% list_var_NPHNLL
    var_NPH_NPHNLL <- list_var_NPH %in% list_var_NPHNLL
    
#    if( any(var_LIN_NPHNLL)){
#      stop("ERROR: some linear variables are also non-linear or non-proportional variables")
#    }

#    if( any(var_NLL_NPHNLL)){
#      stop("ERROR: some non-linear variables are also non-linear-non-proportional variables")
#    }

#    if( any(var_NPH_NPHNLL)){
#      stop("ERROR: some non-proportional variables are also non-linear-non-proportional variable")
#    }
    
                                        # build designs matrices
      des <- ReadDesignFlexrsurv(Terms=Terms, modframe=mf, data=data, rate=rate, Spline_t0=Spline_t0)


# test if all time spline are identical
      if(!is.null(des$X) | !is.null(des$Z)){
        allSpline_T <- c(des$Spline_XT, des$Spline_ZT)
        if(length(allSpline_T) > 1){
          for(sb in allSpline_T){
            if(!identical(sb, allSpline_T[[1]])){
              stop("flexrsurv cannot handle different spline basis for time dependent effects")
            }
          }
        }
        Spline_t <-allSpline_T[[1]] 
      }

                                        # rate
    therate<-des$rate

    
# X0 linear and non linear effects
    X0<-des$X0



# X (NPH effects)
    if(!is.null(des$X)){
      X<-as.matrix(des$X)
      #  Intercept.t of NPH() effets are set in Flersurv:fix.flexrsurv.formula
      #  now build the vector of Intercept_t_NPH for each X var
      # {linear or non linear} and non prop effect are possible; set intercept for NPH effects
      Intercept_t_NPH <- rep(TRUE, length(des$XVars))
      for (i in attr(des$TermsX, "specials")[["NPH"]]){
        thecall <-  match.call(NPH, attr(des$TermsX,"variables")[[i+1]])
        # the variable name is the second argument of the special function
        Intercept_t_NPH[i] <- ifelse( length(thecall[["Intercept.t"]]) == 0 ,
                                     formals(NPH)[["Intercept.t"]],
                                     thecall[["Intercept.t"]] )
      }
    } else {
      X <- NULL
      Intercept_t_NPH <- NULL
    }

# Z (NPHNLL effects


    # get splines for each NPHNLL variable

    Spline_Z <- des$Spline_Z
    if( !is.null(des$Z) ){
      Z<-DesignMatrixNPHNLL(Z=des$Z, listsplinebasis=Spline_Z, timesplinebasis=Spline_t)
    }  else {
      Z <- NULL
    }


  # get initvalues

listinit<- list(gamma0 = init[des$coef2param$gamma0], 
                alpha0 = init[des$coef2param$alpha0],
                beta0  = init[des$coef2param$beta0],
                alpha  = init[des$coef2param$alpha],
                beta   = init[des$coef2param$beta])

# get methods
  
      if( int_meth=="GLM" ){        
        if (is.null(bands)){
          bands <- default_bands(Spline_t)
        }
        method <- list(int_meth=int_meth, bands=bands, optim_meth=optim_meth)
      } else {
        if (is.null(stept)){
          stept <- Max_T /500
        }
        method <- list(int_meth=int_meth, stept=stept, optim_meth=optim_meth)
      }

# fit
    if(ncol(Y)==2){
      # assume Surv(time, event)
        fit<-flexrsurv.ll.fit(X0=X0, X=X, Z=Z, Y=Y,
                              expected_rate=rate,
                              weights=weights,
                              Spline_t0=Spline_t0, Intercept_t0=TRUE,
                              Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
                              bhlink=bhlink,
                              init=listinit,
                              optim.control=optim.control,
                              method=method,
                              vartype = vartype,
                              debug=debug,
                              debug.ll=debug.ll, debug.gr=debug.gr) 




      } else {
      # assume Surv(from, to, event), with possible several observations by subject
        fit<-flexrsurv.ll.fromto.fit(X0=X0, X=X, Z=Z, Y=Y,
                                     expected_rate=rate, 
                                     weights=weights,
                                     Spline_t0=Spline_t0, Intercept_t0=TRUE,
                                     Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
                                     bhlink=bhlink,
                                     init=listinit,
                                     optim.control=optim.control,
                                     method=method,
                                     vartype = vartype,
                                     debug=debug,
                                     debug.ll=debug.ll, debug.gr=debug.gr) 





      }
# returned value  
# buid an object with attributes similares to glm object.
    objfit <- fit


#reorder coefs
    objfit$coefficients <- fit$coefficients[des$param2coef]


    names(objfit$coefficients)[-(1:des$df.T0)] <- des$names_coef_XZ



  
   if( !is.null(dim(fit$var))){
      objfit$var <- (fit$var[des$param2coef, ])[,des$param2coef]  
      dimnames(objfit$var)[[1]] <- names(objfit$coefficients)
      dimnames(objfit$var)[[2]] <- names(objfit$coefficients)
      attr(objfit$var, "type") <- vartype
    }
 if( !is.null(dim(fit$informationMatrix))){
      objfit$informationMatrix <- (fit$informationMatrix[des$param2coef, ])[,des$param2coef]  
      dimnames(objfit$informationMatrix)[[1]] <- names(objfit$coefficients)
      dimnames(objfit$informationMatrix)[[2]] <- names(objfit$coefficients)
      attr(objfit$informationMatrix, "type") <- vartype
    }
  objfit$des <- des
  objfit$terms <- Terms
# don't know if it is necessary?
#  attr(objfit$terms, ".Environment") <- parent.frame()
  objfit$assign <- des$assign
  objfit$assignList <- des$assignList
  
  objfit$na.action <- attr(m, "na.action")

  objfit$optim.control <- optim.control

  objfit$converged <- objfit$conv
  objfit$conv <- NULL

  class(objfit) <- "flexrsurv.mle"

return(objfit)

}
