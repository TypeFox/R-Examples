# latest change done 1-11-13 MS 
# this is an attempt to create a function which 
# fits a GAMLSS distribution using non-linear maximisation like optim()
# in fact here we use MLE() which is a copy of the mle() function of stat4
# The reason for doing this 
#     i) to improve the speed of fitting 
#     ii)  because it would be useful to extent this to Hidden Markov Models
#          and other models as ARMA and GARCH type
# Author Mikis Stasinopoulos 
# TO DO
# i)   weights implementation OK 
# ii)  residuals  OK
# iii) do I need mu.fv and mu.lp?? probably not
# iv)  mu.coef etc OK
# v) if implemented may histDist() should use it (Not working at the moment)
# vi) check all the methods OK
# vii) do I need data argument OK
# viii) vcoc() OK
# needs a summary() function OK
# new TO DO  
# i) If the Hessian is falling may still be able to find parameters (we need to change MLE to allow for that) ok 22-12-11
# ii) we need start.from as argument OK 21-12-11
# latest change 10-7-12
# latest changes are done to correct some of the problems with vcov() function
# only the MLE function is mainly effected
# the optimHess() function is now used for the hessian matrix
# the method optim.proc() in MLE allows to use  optim() or nlminb()

# can we use profile likelihood here?
# not yet but I am working on it now 
#######################################################################################
#names(m1)          
#[5]  "weights"                   
#        "iter"              "method"         
#[13] "converged"       "residuals"       "noObs"          
#[17] "mu.fv"           "mu.lp"            
#[21] "mu.link"              
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# this should be a GAMLSS based mle procedure for fitting distributions
# the main problem is to create the likelihood function automatcally 
# not to be given by the user as in MLE below
# also it would be better to use the link function parameters rather the original
# to avoid boundary problems in the parameters
# this probably can be done by 
#   i)    using the eta.par as parameters in the likelihood fun say LogLikelihood(tparameters)
#   ii)   within the function use the inverse link to go to the original parameters
#   iii)  evalute the likelihood 
#   vi)   after exit transfer back (also use this for the se's)?? 
#---------------------------------------------------------------------------------------- 
#require(gamlss)
gamlssML<-function(y, 
		         family = NO,  
		        weights = NULL, 
		       mu.start = NULL, 
        sigma.start = NULL, 
	         nu.start = NULL, 
		      tau.start = NULL,
	           mu.fix = FALSE,
	        sigma.fix = FALSE,
		         nu.fix = FALSE,
		        tau.fix = FALSE,
               data = NULL, # it need start.from = NULL
         start.from = NULL,
			              ...) 
{
#------------------------------------------------------------------------------------------
# local functions
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
# this is  copy of the stats::mle()
# but it exports a S3 object rather than an S4 object
# see  the examples in mle() for its use
# here we use it as convenient way of calling optim()
# and to cover the case that some parameters can be fixed
#-------------------------------------------------------------------------------------------
# local function taken from the mle() function
# this function is not using hessian at the fitting
#-------------------------------------------------------------------
  MLE <- function (minuslogl, 
                     start = formals(minuslogl), 
                    method = "BFGS", 
  		             fixed = list(),
			    optim.proc = c( "nlminb", "optim"),
			 optim.control = NULL,
			           ...) 
	{
#------------------------------------------------------------------------
	  	  call <- match.call()
	   optim.p <- match.arg(optim.proc)
	  	     n <- names(fixed)
	  fullcoef <- formals(minuslogl)
		if (any(!n %in% names(fullcoef))) 
			stop("some named arguments in 'fixed' are not arguments to the supplied log-likelihood")
 fullcoef[n] <- fixed
		if (!missing(start) && (!is.list(start) || is.null(names(start)))) 
			stop("'start' must be a named list")
   	start[n] <- NULL
		   start <- sapply(start, eval.parent)
		      nm <- names(start)
		      oo <- match(nm, names(fullcoef))
		if (any(is.na(oo))) 
			stop("some named arguments in 'start' are not arguments to the supplied log-likelihood")
		   start <- start[order(oo)]
		      nm <- names(start)
		       f <- function(p) 
		       {
			           l <- as.list(p)
			    names(l) <- nm
			        l[n] <- fixed
			    do.call("minuslogl", l)
		       }
  		 if (length(start))
  		  {
  		  switch(optim.p, 
  		       "nlminb"={
              oout <- nlminb(start = start, objective = f,  control=optim.control)
                   if (oout$convergence > 0) # I took this from Ripley
                   warning("possible convergence problem: optim gave code=", 
                            oout$convergence, " ", oout$message)
       oout$hessian <- optimHess(oout$par, f) #HessianPB(pars=oout$par, fun=f)$Hessian
         value.of.f <- oout$objective
                         },
                 "optim"={
               oout <- optim(par= start, fn = f, method=method, control=optim.control,...)
              #oout <- optim(start, f, method = method, hessian = TRUE, ...)
                   if (oout$convergence > 0) # I took this from Ripley
                   warning("possible convergence problem: optim gave code=", 
                            oout$convergence, " ", oout$message)
        oout$hessian <- optimHess(oout$par, f) #HessianPB(pars=oout$par, fun=f)$Hessian
          value.of.f <- oout$value
                          })    
  		    } 
	   	  else 
			oout <- list(par = numeric(0L), value = f(start))		
		    coef <- oout$par
	  if (length(coef))
  	     {
		  	 vcov <- try(solve(oout$hessian))
            if (any(class(vcov)%in%"try-error")) vcov <- matrix(NA, dim(oout$hessian)[1], dim(oout$hessian)[2])		  	
	     } 				
	else     vcov <- matrix(numeric(0L), 0L, 0L)
     fullcoef[nm] <- coef
           method <- if (optim.p =="nlminb") "nlminb" else method
		#new("mle", call = call, coef = coef, fullcoef = unlist(fullcoef), 
		#    vcov = vcov, min = min, details = oout, minuslogl = minuslogl, 
		#    method = method)
		out <- list(call = call, coef = coef, fullcoef = unlist(fullcoef), 
				vcov = vcov, min = value.of.f, details = oout, minuslogl = minuslogl, 
				method = method)
class(out) <- "MLE"
		out
	}
# end of MLE
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# this is to replicate rqres within gamlss enviroment DS Friday, March 31, 2006 at 10:30
# it is used as in gamlss()
rqres <- function (pfun = "pNO", 
                   type = c("Continuous", "Discrete", "Mixed"),
               censored = NULL,  
                   ymin = NULL, 
                 mass.p = NULL, 
                prob.mp = NULL,
                      y = y,
                         ... )
{ }
body(rqres) <-  eval(quote(body(rqres)), envir = getNamespace("gamlss"))
##---------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# main function starts here
  mlFitcall <- match.call()  #   the function call  
  # if data exit attach them
  #     if (!is.null(data)) {attach(data, name="The_Data_Env"); on.exit(detach("The_Data_Env"))}
      YY <- deparse(substitute(y))     # we need this to allow formulas 
    if (!is.null(data))                # if data exist
    {
     y  <-  if (grepl("~", YY) ) get(as.character(y[2]), envir=as.environment(data)) else      # formula 
                   get(deparse(substitute(y)), envir=as.environment(data))   # non formula
    }
    if (is.null(data))
    {
    y <-  if (grepl("~", YY) ) stop("with formula you need to use the data argument") else y
    }
       fam  <- as.gamlss.family(family)
      fname <- fam$family[[1]] 
       dfun <- paste("d",fname,sep="")
     # pfun <- paste("p",fname,sep="")
        PDF <- eval(parse(text=dfun))
      #	CDF <- eval(parse(text=pfun))
      nopar <- fam$nopar
          N <- length(y)
          w <- if(is.null(weights))    rep(1, N) else weights
            if(any(w < 0)) stop("negative weights not allowed") # 
## extracting now the y and the binomial denominator in case we use BI or BB
        if(any(fam$family%in%.gamlss.bi.list)) 
           { 
            if (NCOL(y) == 1) # binary
               {
                   y <- if (is.factor(y))  y != levels(y)[1] else y
                  bd <- rep(1, N)
                  if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
               } 
            else if (NCOL(y) == 2) 
               {
                 if (any(abs(y - round(y)) > 0.001)) 
                   {
                    warning("non-integer counts in a binomial GAMLSS!")
                    }
                 bd <- y[,1] + y[,2]
                  y <-  y[,1]
              if (any(y < 0 | y > bd)) stop("y values must be 0 <= y <= N") # MS Monday, October 17, 2005 
               } 
              else stop(paste("For the binomial family, y must be", 
                  "a vector of 0 and 1's or a 2 column", "matrix where col 1 is no. successes", 
                  "and col 2 is no. failures"))
            }
# multinomial checking
            else if(any(fam$family%in%.gamlss.multin.list))
                {
                  y <- if(is.factor(y))   unclass(y)
                       else y
               } 
# here is the place to  start from ACTION HERE
# check whether start.from is null
if (!is.null(start.from)) # MS 21-12-11
 {
       if (is.gamlss(start.from)) # if not check whether model or vector
       {
       	if ("mu"%in%names(fam$parameters))    {mu.start <- fitted(start.from, "mu")[1]}
       	if ("sigma"%in%names(fam$parameters)) {sigma.start <- fitted(start.from, "sigma")[1]}
       	if ("nu"%in%names(fam$parameters))    {nu.start <- fitted(start.from, "nu")[1]}
       	if ("tau"%in%names(fam$parameters))   {tau.start <- fitted(start.from, "tau")[1]}
       }
       else 
       {
       	if (is.numeric(start.from))
       	 {
       	if (nopar<length(start.from)) stop("start.fom need to be as big as the number of parameters in the distribution")
       	if ("mu"%in%names(fam$parameters))    {mu.start <- start.from[1]}
       	if ("sigma"%in%names(fam$parameters)) {sigma.start <- start.from[2]}
        if ("nu"%in%names(fam$parameters))    {nu.start <- start.from[3]}
        if ("tau"%in%names(fam$parameters))   {tau.start <- start.from[4]}
       	 }
       	 else stop("start.from is a non numeric variable")
       }
 }
##  get the initial values if are not set by the user
           if ("mu"%in%names(fam$parameters))
           {
                  mu <- if(is.null(mu.start))  mean(eval(fam$mu.initial, list(y=y)))
		       	       else mu.start[1] 
              eta.mu <- fam$mu.linkfun(mu)
			#eta.par <- c(eta.mu)          
           }
           if ("sigma"%in%names(fam$parameters))
           {
               sigma <- if(is.null(sigma.start)) mean(eval(fam$sigma.initial, list(y=y, mu=mu)))
			            else sigma.start[1]
           eta.sigma <- fam$sigma.linkfun(sigma)
			#eta.par <- c(eta.mu, eta.sigma)
           }
           if ("nu"%in%names(fam$parameters))
           {
			      nu <- if(is.null(nu.start))  mean(eval(fam$nu.initial, list(y=y, mu=mu, sigma=sigma)))
			            else nu.start[1]
              eta.nu <- fam$nu.linkfun(nu)
			#eta.par <- c(eta.mu, eta.sigma, eta.nu)
           }
           if ("tau"%in%names(fam$parameters))
           {  
			     tau <- if(is.null(tau.start)) mean(eval(fam$tau.initial, list(y=y, mu=mu, sigma=sigma, nu=nu)))
			           else tau.start[1]
           	 eta.tau <- fam$tau.linkfun(tau)
			#eta.par <- c(eta.mu, eta.sigma, eta.nu, eta.tau)
           }
## whether to fix parameters
   fixed <- list()
   if (mu.fix)      fixed <- c(fixed, eta.mu=mu[1])
   if (sigma.fix)   fixed <- c(fixed, eta.sigma=sigma[1])
   if (nu.fix)      fixed <- c(fixed, eta.nu=nu[1])
   if (tau.fix)     fixed <- c(fixed, eta.tau=tau[1])
                  noFixed <- sum(c(mu.fix, sigma.fix, nu.fix, tau.fix))
## define the likelihood and find the maximum
    switch(nopar, 
			{# one parameter 
				ll1 <- function(eta.mu)
				{
					mu <- fam$mu.linkinv(eta.mu)
					if(any(fam$family%in%.gamlss.bi.list)) -sum(w*PDF(y, mu=mu, bd=bd, log=TRUE)) # BI
          else -sum(w*PDF(y, mu=mu, log=TRUE))# other
				}
			   fit <-	MLE(ll1, start=list(eta.mu=eta.mu), fixed=fixed, ...)
			},
			{# two paremeters
				ll2 <- function(eta.mu, eta.sigma)
				{
					mu <- fam$mu.linkinv(eta.mu)
				 sigma <- fam$sigma.linkinv(eta.sigma)
					if(any(fam$family%in%.gamlss.bi.list))  -sum(w*PDF(y, bd=bd, mu=mu, sigma=sigma, log=TRUE))  else   
					-sum(w*PDF(y, mu=mu, sigma=sigma, log=TRUE))
				}
				fit <-	MLE(ll2, start=list(eta.mu=eta.mu, eta.sigma=eta.sigma), fixed=fixed, ...)
			},
			{# three parameters
				ll3 <- function(eta.mu, eta.sigma, eta.nu)
				{
					   mu <- fam$mu.linkinv(eta.mu)
					sigma <- fam$sigma.linkinv(eta.sigma)
					   nu <- fam$nu.linkinv(eta.nu)
			if(any(fam$family%in%.gamlss.bi.list)) -sum(w*PDF(y, bd=bd, mu=mu, sigma=sigma, nu=nu, log=TRUE)) else
               -sum(w*PDF(y, mu=mu, sigma=sigma, nu=nu, log=TRUE))
				}
				fit <-	MLE(ll3, start=list(eta.mu=eta.mu, eta.sigma=eta.sigma, eta.nu=eta.nu), fixed=fixed, ...)
				
			},
			{# four parameters
				ll4 <- function(eta.mu, eta.sigma, eta.nu, eta.tau)
				{
				    	mu <- fam$mu.linkinv(eta.mu)
				     sigma <- fam$sigma.linkinv(eta.sigma)
					    nu <- fam$nu.linkinv(eta.nu)
                       tau <- fam$tau.linkinv(eta.tau)
					-sum(w*PDF(y, mu=mu, sigma=sigma, nu=nu, tau=tau, log=TRUE))
				}
				fit <-	MLE(ll4, start=list(eta.mu=eta.mu, eta.sigma=eta.sigma, eta.nu=eta.nu, eta.tau=eta.tau), fixed=fixed, ...)
			}
	  ) # end of switch
# saving things
  df.fit <- nopar-noFixed
     out <-list(family = fam$family,  
			parameters =  as.character(names(fam$par)), 
			      type = fam$type, 
				  call = mlFitcall, 
				     y = y,
			   weights = w,
	        G.deviance = 2*fit$min, 
			         N = N,  
				df.fit = df.fit, 
		   df.residual = N-df.fit, 
				   aic = 2*fit$min+2*nopar,
				   sbc = 2*fit$min+log(N)*nopar,
				method = fit$method,
		   		  vcov = fit$vcov, 
		        Allpar = fit$coef)
    if ("mu"%in%names(fam$parameters))
		   {
			            mu <- if (mu.fix)  fam$mu.linkinv(eta.mu) else fam$mu.linkinv(fit$coef["eta.mu"])
			     names(mu) <- "mu"
		   mu.coefficients <- fit$coef["eta.mu"]
    names(mu.coefficients) <- "mu.coefficients"
				       out <- c(out, mu, mu.coefficients )
			   out$mu.link <- fam$mu.link
		   }
		   ## Output for sigma model: ---------------------------------------------------------------
		   if ("sigma"%in%names(fam$parameters))
		   {
			        sigma <- if (sigma.fix)  fam$mu.linkinv(eta.sigma) else fam$sigma.linkinv(fit$coef["eta.sigma"])
		     names(sigma) <- "sigma"
 	   sigma.coefficients <- fit$coef["eta.sigma"]
names(sigma.coefficients) <- "sigma.coefficients"	
				      out <- c(out, sigma, sigma.coefficients )
		   out$sigma.link <- fam$sigma.link
		   }
		   ##  output for nu ------------------------------------------------------------------------
		   if ("nu"%in%names(fam$parameters))
		   {
			     nu <- if (nu.fix)  fam$mu.linkinv(eta.nu) else fam$nu.linkinv(fit$coef["eta.nu"])
				names(nu) <- "nu"
				nu.coefficients = fit$coef["eta.nu"]
				names(nu.coefficients) <- "nu.coefficients"			    
				 out <- c(out, nu, nu.coefficients )
				 out$nu.link <- fam$nu.link

		   }
		   ##  output for tau -----------------------------------------------------------------------
		   if ("tau"%in%names(fam$parameters))
		   {
			   tau <- if (tau.fix)  fam$mu.linkinv(eta.tau) else fam$tau.linkinv(fit$coef["eta.tau"])
				names(tau) <- "tau"
				tau.coefficients = fit$coef["eta.tau"]
				names(tau.coefficients) <- "tau.coefficients"			    
				 out <- c(out, tau, tau.coefficients )
				 out$tau.link <- fam$tau.link
		   }
		    out$residuals <- eval(fam$rqres)
		        out$rqres <- fam$rqres
		  if (!is.null(data) ) out$call$data <- substitute(data)
     class(out) <- c("gamlssML", "gamlss")
	 out
 }                                        
######################################################################################
# methods for gamlssML    
######################################################################################
fitted.gamlssML<-function (object, what = c("mu", "sigma", "nu", "tau"), parameter= NULL, ... ) 
{
what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
if (! what%in%object$par) stop(paste(what,"is not a parameter in the gamlss object","\n"))
x <- rep(object[[what]], object$N)
x
}
######################################################################################
vcov.gamlssML<-function (object, type = c("vcov", "cor", "se",  "all"), ... ) 
{
type <- match.arg(type)
switch(type,
       "se"={x <- sqrt(diag(object$vcov))},
       "cor"={x <- cov2cor(object$vcov)},
       "vcov"={x<- object$vcov},
       	"all"={x <-list(vcov=object$vcov, cor=cov2cor(object$vcov), se=sqrt(diag(object$vcov)))}
       	)
x
}
######################################################################################
# new at the 6-8-11 DS  
summary.gamlssML  <- function (object, digits = max(3, getOption("digits") - 3), ...) 
{   
	
     digits <- max(3, getOption("digits") - 3)
      cat("*******************************************************************")
    cat("\nFamily: ", deparse(object$family), "\n") 
    cat("\nCall: ",  deparse(object$call),  "\n", fill=TRUE)
    cat("Fitting method:", deparse(object$method), "\n\n") 
    est.disp <- FALSE
       # df.r <- object$noObs - object$mu.df
            coef  <- object$Allpar
          se.coef <- vcov(object, "se")
             tval <- coef/se.coef
          matcoef <- cbind(coef, se.coef, tval, 2*(1-pnorm(abs(tval))))
dimnames(matcoef) <- list(names(tval), c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
cat("\nCoefficient(s):\n")
printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
     cat("\n Degrees of Freedom for the fit:", object$df.fit, "Residual Deg. of Freedom  ", 
        object$df.residual, "\n")
    cat("Global Deviance:    ", format(signif(object$G.deviance)), 
        "\n            AIC:    ", format(signif(object$aic)), "\n            SBC:    ", 
        format(signif(object$sbc)), "\n")
    invisible(object)
}
# methods are finish here  
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------


   
 
 


