# for Vlasi look for ##  comments for the family ##
# last change 2-4-11
# garmaFit() is a function based on the Garch11fit() of the fGarch pacakge 
# it uses a non linear maximation routine to fit a GARMA(1,1) model
# it is a usefull application of the recursive use of filter() 
# TO DO
# 
# 1) Can generalised for more general ARMA model arma(p,q)? yes  OK
# 2) for more general distributions yes OK
# 3) needs se error for the initial parametets so we can get lower and upper bounds: maybe we dont need this set to -Inf and Inf
# 4) residuals() OK
# 5) test different link functions for g2(y*) 
# 6) fix to allows 3 and 4 parameter distributions
# 7) fitted()   OK
# 8) deviance() OK
# 9) AIC()      OK
#10) vcov()     OK
#11) summary()  OK
#12) plot()     OK
#13) print()    ok
#14) coef()     ok
#15) wp()
# 16) fixed parameters
#source("/Users/stasinom/Documents/gamlss/projects/TIME-SERIES/createLags.R")
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#require(gamlss)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
garmaFit <- function(formula = formula(data), 
		                   order = c(0,0),
		                 weights = NULL,
	                      data = sys.parent(),
	                    family = NO(),  ##  declaring the family ##
		                   alpha = 0.1,
				           phi.start = NULL,
				         theta.start = NULL, 
			                  tail = max(order),
		                 control = list())
{
#-------------------------------------------------------------------------------
# local functions
# i)   rqres
# ii)  hessian not needed anymore 
# iii) createLags
# iv)  garmaLL
#-------------------------------------------------------------------------------
# (i)
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
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
# the hessian
#  (ii)
#-------------------------------------------------------------------------------
	hessian<-function(par)
	{
		npar <- length(par)
		epsilon <- 0.0001 * par
		Hessian = matrix(0, ncol = npar, nrow = npar)
		for (i in 1:npar) 
		{
			for (j in 1:npar) 
			{
				x1. <- x2. <- x3. <- x4. <- par
				x1.[i] <- x1.[i] + epsilon[i]; x1.[j] <- x1.[j] + epsilon[j] 
				x2.[i] <- x2.[i] + epsilon[i]; x2.[j] <- x2.[j] - epsilon[j]
				x3.[i] <- x3.[i] - epsilon[i]; x3.[j] <- x3.[j] + epsilon[j]
				x4.[i] <- x4.[i] - epsilon[i]; x4.[j] <- x4.[j] - epsilon[j]
				Hessian[i, j] <- (garmaLL(x1.)-garmaLL(x2.)-garmaLL(x3.)+garmaLL(x4.))/(4*epsilon[i]*epsilon[j])
			}
		}
		Hessian
	}
##------------------------------------------------------------------------------
# (iii)
##------------------------------------------------------------------------------
createLags <- function(y,lag=1, omit.na=FALSE)
{
  ## local function   
  lag1 <- function(y)
  {
    l <- c(NA, y[1:(length(y)-1)])
    l
  }
  ## main function starts here
  d <- matrix(0, nrow=length(y), ncol=lag)
  d[,1] <- lag1(y)
  yname <- deparse(substitute(y))
  cname <- paste("lag1", yname, sep="")
  if (lag==1) 
  {names(d[,1]) <- c("lag1")
   if (omit.na) d <-na.omit(d)
   return(d)
  }
  for(i in 2:lag) 
  {
    d[,i] <- lag1(d[,i - 1])
    cname <- c(cname, paste("lag", paste(yname,i, sep=""), sep=""))
  }
  colnames(d ) <- cname
  if (omit.na) d <-na.omit(d)
  d
}
#-------------------------------------------------------------------------------
# (iv)
##------------------------------------------------------------------------------
	garmaLL <- function(parm,  save=FALSE) 
	{	
		           beta <- parm[1:l.beta]	
		             Xb <- X%*%beta
		switch(case, 
				{ # case 1
			      start  <- l.beta + 1
				    finish <- start + (order[1]-1)
				       phi <-  parm[start:finish]
				  finish   <- finish + 1
				    g2ylag <- createLags(g2y, order[1])
				     Xblag <- createLags(Xb, order[1])
				    g2y_Xb <- g2ylag-Xblag
				       eta <- Xb +  g2y_Xb%*%phi
				        mu <- ifelse(is.na(eta), alpha, fam$mu.linkinv(eta)) 
				},
				{ # case 2
					
				  	start  <- l.beta + 1
		  		 	finish <- start+ (order[2]-1)
			  		 theta <- parm[start:finish]
				    finish <- finish + 1
				  	g2ylag <- createLags(g2y, order[2])     
					    some <- Xb+g2ylag%*%theta
					    some <- ifelse(is.na(some),0,some) 
					     eta <- filter(some, -theta, "r", init =rep(0,order[2]))
					    mu <- ifelse(is.na(eta), alpha, fam$mu.linkinv(eta)) 	
				},
				{ # case 3
					#print(parm)
					#if(any(is.na(parm))) browser()
				  	  start <- l.beta + 1
			  		finish1 <- start + (order[1]-1)
			  		finish2 <- finish1+1+ (order[2]-1)
			  		   phi  <- parm[start:finish1]
				  	  theta <- parm[(finish1+1):finish2]
			  		 finish <- finish2+1    
			 	   g2ylagar <- createLags(g2y, order[1])    
				   g2ylagma <- createLags(g2y, order[2])
				    Xblagar <- createLags(Xb, order[1])
				    Xblagma <- createLags(Xb, order[2])
				  	 g2y_Xb <-  g2ylagar-Xblagar
				  	#print(phi)
				  	# print(theta)
				  	 #print(head(cbind( g2y_Xb, g2ylagar, g2ylagma, Xblagar, Xblagma)))
					     some <- Xb+g2y_Xb%*%phi+g2ylagma%*%theta
					     some <- ifelse(is.na(some),0,some) 
					  # print(head(some))					  
					     eta <- filter(some, -theta, "r", init =rep(0, order[2])) 
			#tryCatch(eta <- filter(some, -theta, "r", init =rep(0, order[2])), finally=print(head(eta)))
					     mu <- ifelse(is.na(eta), alpha, fam$mu.linkinv(eta)) 		  
				}
		       )
			   switch(nopar, ##  this is where swiches between different family with different distributions ##
					   {# one parameter 
						   llh <- if (BItrue) -sum(w*PDF(y, mu=mu, bd=bd, log=TRUE), na.rm = TRUE)
                      else -sum(w*PDF(y, mu=mu,  log=TRUE), na.rm = TRUE)
					   }, 
					   {# two  parameters 
						 sigma <- parm[finish] 
						   llh <-  -sum(w*PDF(y, mu=mu, sigma=sigma, log=TRUE), na.rm = TRUE)
					   },
					   {# three parameters
						 sigma <- parm[finish] 
						    nu <- parm[finish+1]
						   llh <- -sum(w*PDF(y, mu=mu, sigma=sigma, nu=nu, log=TRUE), na.rm = TRUE)
						   
					   },
					   {# four parameters
						 sigma <- parm[finish] 
						    nu <- parm[finish+1]
						   tau <-parm[finish+2]
						   llh <-  -sum(w*PDF(y, mu=mu, sigma=sigma, nu=nu, tau=tau, log=TRUE),  na.rm = TRUE)
             
					   }
			         )
		if (save) return(list(lik=llh, mu=mu))
	    #if(is.na(llh)) browser()
		llh 
	}	
#-------------------------------------------------------------------------------
# the main function starts here
#-------------------------------------------------------------------------------
	garmacall <- match.call()  #   the function call
    if (order[1]<=0&&order[2]<=0)  case <- 0 
    if (order[1]> 0&&order[2]<=0)  case <- 1 
	  if (order[1]<=0&&order[2] >0)  case <- 2	
	  if (order[1]> 0&&order[2] >0)  case <- 3
## starting values for beta
## fitting a gamlss model
## possibly with weights ????? YES but we need length(y) here 
     m0 <- gamlss(formula, family=family, data=data, trace=FALSE) #
     cat("deviance of linear model= ", deviance(m0),"\n")
     if (case==0) return(m0) # stop if case 0
# get the initial betas
       beta <- coef(m0)
   # here we need the se of  coef so we can create lower and upper bounds???
   #  se.coef <- c(0.09662, 0.13091, 0.13968, 0.13302, 0.13437) 
   #  se.coef <- vcov(m0, "se") NEEDS ATTENTION
	   l.beta <- length(beta)
	 #  beta.se <- [1:l.beta]
          y <- m0$y # n
	        X <- m0$mu.x # the X matrix
	      fam <- as.gamlss.family(family)
      nopar <- fam$nopar
          N <- length(y)
## extracting now the y and the binomial denominator in case we use BI or BB
if(any(fam$family%in%gamlss::.gamlss.bi.list)) 
{ BItrue <- TRUE
      bd <- m0$bd 
}  
#   if (NCOL(y) == 1) # binary
#   {
#     y <- if (is.factor(y))  y != levels(y)[1] else y
#     bd <- rep(1, N)
#     if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
#   } 
#   else if (NCOL(y) == 2) 
#   {
#     if (any(abs(y - round(y)) > 0.001)) 
#     {
#       warning("non-integer counts in a binomial GAMLSS!")
#     }
#     bd <- y[,1] + y[,2]
#     y <-  y[,1]
#     if (any(y < 0 | y > bd)) stop("y values must be 0 <= y <= N") # MS Monday, October 17, 2005 
#   } 
#   else stop(paste("For the binomial family, y must be", 
#                   "a vector of 0 and 1's or a 2 column", "matrix where col 1 is no. successes", 
#                   "and col 2 is no. failures"))
# }
	#if (!is.null(data)) {attach(data); on.exit(detach(data))}
   ##  the next 6 lines are needed for the  family ##
      fname <- fam$family[[1]] 
       dfun <- paste("d",fname,sep="")
       pfun <- paste("p",fname,sep="")
	      PDF <- eval(parse(text=dfun))
	      CDF <- eval(parse(text=pfun))
	# depending on the link function
	# I am not sure I cover all possibilities what about own???
	# also the logit needs testing ?????
	     if (fam$mu.link=="identity") ystar <- y 
	     if (fam$mu.link=="log")      ystar <- pmax(y,alpha)
	     if (fam$mu.link=="logit")    ystar <- pmin(pmax(y, alpha), bd-alpha)/bd# ifelse(y==0, alpha, ifelse( y==1, 1-alpha, y) )
	     if (fam$mu.link=="inverse")  ystar <- y
	     g2y <- fam$mu.linkfun(ystar) 
	# I will need here ar and ma 
    tailoff <- max(order)
	    mtail <- max(tailoff, tail)
	        w <- if(is.null(weights))    rep(1, N) else weights
	          if(any(w < 0)) stop("negative weights not allowed") # 
 w[1:mtail] <- 0	# weight out the tail  
# get the initial values if are not set by the user
##  you will need  something like this for the general family probably you need for nu and tau only##
	if ("mu"%in%names(fam$parameters))
	{
    		params <- c(beta = beta) 
   lowerBounds <- c(beta = rep(-Inf, l.beta)) # those bounds need checking ???
   upperBounds <- c(beta = rep( Inf, l.beta))
		 switch(case, 
				 {
				 	
				 	if (length(phi.start)!=order[1]&&!is.null(phi.start)) 
				 	        stop("phi.start should have ", order[1], " elements \n" )
				      phi <- if (is.null(phi.start)) runif(order[1]) else phi.start           
			       params <- c(params,       phi = phi)
			  lowerBounds <- c(lowerBounds,  phi = rep(-1,order[1]))
			  upperBounds <- c(upperBounds,  phi = rep(1, order[1])) 
			      },
		          {
		          	if (length(theta.start)!=order[2]&&!is.null(theta.start)) 
				 	        stop("theta.start should have ", order[2], " elements \n" )
				    theta <- if (is.null(theta.start)) runif(order[2]) else theta.start  
			       params <- c(params,      theta= theta)
			  lowerBounds <- c(lowerBounds, theta = rep(-1,order[2]))
		      upperBounds <- c(upperBounds, theta = rep(1, order[2])) 	  
			      },
				  {
				  	if (length(phi.start)!=order[1]&&!is.null(phi.start)) 
				 	        stop("phi.start should have ", order[1], " elements \n" )
				 	if (length(theta.start)!=order[2]&&!is.null(theta.start)) 
				 	        stop("theta.start should have ", order[2], " elements \n" )
				 	   phi <- if (is.null(phi.start)) runif(order[1]) else phi.start 
				 	 theta <- if (is.null(theta.start)) runif(order[2]) else theta.start           
				   params <- c(params,       phi = phi,   theta = theta)
		      lowerBounds <- c(lowerBounds,  phi = rep(-1,order[1]),  theta = rep(-1,order[2]))
		      upperBounds <- c(upperBounds,  phi = rep(1, order[1]),  theta = rep(1,order[2])) 	  
				  }
		       )
    }
	if ("sigma"%in%names(fam$parameters))
	{
		     sigma <- fitted(m0, "sigma")[1] 
	  names(sigma) <- ""
		#sigma.coef <- coef(m0, "sigma")
		#  sigma.se <- se.coef[l.beta+1]
		# This needs fixxing  
		#lowersigma <-  fam$smgma.linkinv()
		params <- c(params, sigma = sigma )
		lowerBounds <- c(lowerBounds, sigma =  0.001)
		upperBounds <- c(upperBounds, sigma =  1e10) 
	}
	if ("nu"%in%names(fam$parameters))  ##  probably ypu need something similar ##
	{
		nu <- fitted(m0, "nu")[1] 
		names(nu) <- ""
		#sigma.coef <- coef(m0, "sigma")
		#  sigma.se <- se.coef[l.beta+1]
		# This needs fixxing  
		#lowersigma <-  fam$smgma.linkinv()
		params <- c(params, nu = nu )
		lowerBounds <- c(lowerBounds, nu = -Inf)# I have to think about this ?
		upperBounds <- c(upperBounds, nu =  Inf) 
	}
	if ("tau"%in%names(fam$parameters)) ##  probably ypu need something similar ##
	{  
		tau <- fitted(m0, "tau")[1] 
		names(tau) <- ""
		params <- c(params, tau = tau )
		lowerBounds <- c(lowerBounds, tau = -Inf)# I have to think about this ?
		upperBounds <- c(upperBounds, tau =  Inf)
	} 
	     # S <- 1e-6
	
#  Estimate Parameters 
##  now we have to make sure that the likelihood can cope with different famillies so see garmaLL ##
#	fit <- nlminb(start = params, objective = garmaLL, lower = lowerBounds,  upper = upperBounds, control = control)
#  browser()
# optim(par = params, fn = garmaLL, lower = lowerBounds,  upper = upperBounds, control = control, method="L-BFGS-B")  
 				   fit <- try( nlminb(start = params, objective = garmaLL, lower = lowerBounds,  upper = upperBounds, control = control))
					   if (any(class(fit)%in%"try-error"))
                 { 
                 cat("OOPS IT FAILED: let us try once more \n")	
                  fit <- do.call("garmaFit", args=as.list(garmacall[-1]))
                  return(fit)
                 } 
   
    
    # fit <-  optim(params, garmaLL, method = c("L-BFGS-B"),lower = lowerBounds, upper = upperBounds, control = list(), hessian = FALSE)
   fit1 <- garmaLL(fit$par, save=TRUE)  
 df.fit <- length(fit$par)
Hessian <- hessian(fit$par)   
# Step 6: Create and Print Summary Report:
#          se.coef <- sqrt(diag(solve(Hessian)))
#             tval <- fit$par/se.coef
#          matcoef <- cbind(fit$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
#dimnames(matcoef) <- list(names(tval), c(" Estimate",
#" Std. Error", " t value", "Pr(>|t|)"))
#cat("\nCoefficient(s):\n")
#printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
# preparing for output
cat("deviance of  garma model= ", 2*fit$objective,"\n")
if (2*fit$objective>deviance(m0)) warning("There is a problem with the fitted GARMA model here \n")
out <-list(family = fam$family,  
       parameters = as.character(names(fam$par)), 
		     type = fam$type, 
		     call = garmacall, 
		        y = y,
	      weights = w,
       G.deviance = 2*fit$objective, 
	       	    N = N,  
	       df.fit = length(fit$par), 
      df.residual = N-df.fit, 
	 	      aic = 2*fit$objective+2*df.fit,
		      sbc = 2*fit$objective+log(N)*df.fit,
	       method = "nlminb",
		     vcov = solve(Hessian),
		     coef = fit$par) 
if (   "mu"%in%names(fam$parameters))      	
{		
     		      mu <- as.vector(fit1$mu)
	    # mu[1:mtail] <- NA
 		      out$mu <- mu
			  if (fam$nopar>1) 
		       {
				   final <-which(names(fit$par)=="sigma")-1
     out$mu.coefficients <- fit$par[1:final]
		       }
		      else 
			  {  
				  out$mu.coefficients <- fit$par
			  }
	  }
if ("sigma"%in%names(fam$parameters)) 
      {
	  	        sigma <- fit$par[ "sigma"]
	        out$sigma <- sigma
out$sigma.coefficients<- fam$sigma.linkinv(sigma)
	  }
if (   "nu"%in%names(fam$parameters))    
	  {
		         nu <- fit$par[ "nu"]
	         out$nu <- nu
out$nu.coefficients <- fam$nu.linkinv(nu)		 
	  }
if (  "tau"%in%names(fam$parameters))
      {
	             tau <- fit$par[ "tau"]
	         out$tau <- tau
out$tau.coefficients <- fam$tau.linkinv(tau)	
      }
	out$residuals <- eval(fam$rqres)
	    out$rqres <- fam$rqres
	class(out) <-c("garma", "gamlss")
	out
}
######################################################################################
#-------------------------------------------------------------------------------------
# garmafit finish here
#-------------------------------------------------------------------------------------
######################################################################################
vcov.garma<-function (object, type = c("vcov", "cor", "se",  "all"), ... ) 
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
######################################################################################
fitted.garma<-function (object, what = c("mu", "sigma", "nu", "tau"), ... ) 
{
	what <- match.arg(what)
	if (! what%in%object$par) stop(paste(what,"is not a parameter in the gamlss object","\n"))
	x <- if (what=="mu") as.vector(object$mu)
	     else rep(object[[what]], object$N)
	x
}
######################################################################################
######################################################################################
print.garma <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nFamily: ", deparse(x$family), "\n")
    cat("Fitting method:", deparse(x$method), "\n")
    cat("\nCall: ", deparse(x$call), "\n", fill = TRUE)
    cat("Mu Coefficients")
    if (is.character(co <- x$contrasts)) 
        cat("  [contrasts: ", apply(cbind(names(co), co), 1, 
            paste, collapse = "="), "]")
    cat(":\n")
    if ("mu" %in% x$parameters) {
        print.default(format(coef(x, "mu"), digits = digits), 
            print.gap = 2, quote = FALSE)
    }
    if ("sigma" %in% x$parameters) {
        cat("Sigma Coefficients:\n")
        print.default(format(coef(x, "sigma"), digits = digits), 
            print.gap = 2, quote = FALSE)
    }
    if ("nu" %in% x$parameters) {
        cat("Nu Coefficients:\n")
        print.default(format(coef(x, "nu"), digits = digits), 
            print.gap = 2, quote = FALSE)
    }
    if ("tau" %in% x$parameters) {
        cat("Tau Coefficients:\n")
        print.default(format(coef(x, "tau"), digits = digits), 
            print.gap = 2, quote = FALSE)
    }
    cat("\n Degrees of Freedom for the fit:", x$df.fit, "Residual Deg. of Freedom  ", 
        x$df.residual, "\n")
    cat("Global Deviance:    ", format(signif(x$G.deviance)), 
        "\n            AIC:    ", format(signif(x$aic)), "\n            SBC:    ", 
        format(signif(x$sbc)), "\n")
    invisible(x)
}
######################################################################################
######################################################################################
coef.garma <- function (object, what = c("mu", "sigma", "nu", "tau"), ...) 
{
    what <- match.arg(what)
    if (!what %in% object$par) 
        stop(paste(what, "is not a parameter in the object", 
            "\n"))
    x <- object[[paste(what, "coefficients", sep = ".")]]
    x
}#
#####################################################################################
#####################################################################################
summary.garma <- function (object, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nFamily: ", deparse(object$family), "\n")
    cat("Fitting method:", deparse(object$method), "\n")
    cat("\nCall: ", deparse(object$call), "\n", fill = TRUE)
   # cat("Mu Coefficients")
   # if (is.character(co <- object$contrasts)) 
   # cat("  [contrasts: ", apply(cbind(names(co), co), 1, 
   #         paste, collapse = "="), "]")  
   # Step 6: Create and Print Summary Report:
          se.coef <- sqrt(diag(object$vcov))
             tval <- object$coef/se.coef
          matcoef <- cbind(object$coef, se.coef, tval, 2*(1-pnorm(abs(tval))))
dimnames(matcoef) <- list(names(tval), c(" Estimate",
" Std. Error", " t value", "Pr(>|t|)"))
cat("\nCoefficient(s):\n")
printCoefmat(matcoef, digits = 6, signif.stars = TRUE)

    cat("\n Degrees of Freedom for the fit:", object$df.fit, "Residual Deg. of Freedom  ", 
        object$df.residual, "\n")
    cat("Global Deviance:    ", format(signif(object$G.deviance)), 
        "\n            AIC:    ", format(signif(object$aic)), "\n            SBC:    ", 
        format(signif(object$sbc)), "\n")
    invisible(object)
}
#####################################################################################
#####################################################################################
