# $Id: sfa.R 156 2015-07-08 13:34:15Z b002961 $
# \encoding{latin1}


sfa <- function(x, y, beta0=NULL, lambda0=1, resfun=ebeta, 
                TRANSPOSE = FALSE, DEBUG=FALSE, 
                control = list(), hessian=2)  
{

# Funktion: beregner minus loglikelihood
loglik <- function(parm) {
   N <- dim(x)[1] 
   K <- dim(x)[2] + 1
   beta <- parm[1:K]
   lambda <- parm[K+1]
   # s <- parm[K+2]
   e <- resfun(x,y,beta)   
   s <- sum(e^2)/length(e)
   z <- -lambda*e/sqrt(s)
   pz = pmax(pnorm(z),1e-323) # undgaa at der skal tages log af 0, ssh
                              # er altid positiv fordi stoetten er hele 
                              # aksen, afrunding kan saette den til rent 0
   l <- N/2*log(pi/2) +N/2*log(s) -sum(log(pz)) +N/2.0 
   # cat(parm,":  ",l,"\n")
   # print(l)
   return(l)
} # loglik

   if ( class(x) == "data.frame" )  {
      x <- as.matrix(x)
   }

   if ( !is.matrix(x) )  {
      print("Input 'x' must be a matrix",quote=F)
      return(print("Function 'sfa' stops",quote=F))
   }
   if ( is.matrix(y) && ncol(y) > 1 )  {
      stop("Only one column (one variale) allowed for argument y in sfa")
   }
   # Fjern observationer hvor der er NA'er
   cpc <- complete.cases(cbind(x,y))
   if ( sum(!cpc) > 0 )  {
      x <- x[cpc,,drop=FALSE]
      y <- y[cpc]
   }

   if ( TRANSPOSE ) {
      x <- t(x)
      y <- t(y)
   }
   K <- dim(x)[2] + 1
   if ( nrow(x) != length(y) )  {
      stop("Number of observations on x and y differ")
   }

if ( missing(beta0) )  {
   # 1. OLS estimate
   m <- lm( y ~ x)      # OLS estimate; linear model
	# print(logLik(m))
   beta0 <- m$coef
   # sigma0 <- deviance(m)/dim(x)[1]   # estimate of variance, ML estimate
	# print(dim(x))
	# print(loglik(c(m$coef,deviance(m)/dim(x)[1],0)))
}


   # 2. Minimization of minus log likelihood
   parm = c(beta0,lambda0)
   # print(loglik(parm))
   #print("For control er sat"); print(control)
   if ( missing(control) )  {
		#print("control er missing")
		#control["maxeval"] <- 1000
		control["stepmax"] <- .1
   } else {
   	# Hvis 'control' er sat skal det sikres at stepmax har en lille vaerdi
   	# hvis den ikke er sat i 'control'; er den sat i 'control' bliver den vaerdi
   	# benyttet paa brugerens egen risiko.
		#print("control har vaerdier")
		#print(control)
		kontrol <- list()
		kontrol["stepmax"] <- .1
		for ( i in names(control) )  {
			kontrol[i] <- control[i]
		}
		control <- kontrol
   }
   if (DEBUG) {print("Efter control er sat"); print(control)}
   if (DEBUG) print("ucminf bliver kaldt")
   o <- ucminf(parm, loglik, control=control, hessian=hessian)
   if (DEBUG) {
      print("ucminf er slut")
      print(paste("Antal funktionskald",o$info["neval"]))
      print(o$info)
      if ( hessian!= 0 & !is.null(o$hessian) )  {
      	print("Hessian:")
      	print(o$hessian)
      }
      if ( hessian!= 0 & !is.null(o$invhessian) )  {
      	print("Inverse hessian")
   		print(o$invhessian)
   	}
   }

   if ( o$convergence < 0 )  {
      warning("Not necessarily converged, $convergence = ",
                   o$convergence,"\n",o$message)
   }
   if ( DEBUG & o$convergence > 0 )  {
      warning("Converged, $convergence = ", o$convergence,"\n",o$message)
   }
   if (o$par[K+1] < 0)  {
      warning("lambda is negative.\n",
         "  This could indicate that there is no inefficiency,\n",
         "  or that the model is misspecified." )
   }
   sf <- o
   class(sf) <- "sfa"
   K <- dim(x)[2] + 1
   sf$beta <- o$par[1:K]
   sf$coef <- o$par[1:K] 
   names(sf$coef) <- names(o$par[1:K])
   e <- resfun(x,y,sf$beta)   
   sf$residuals <- e
   sf$fitted.values <- y -e
   sf$lambda <- o$par[K+1]
   sf$sigma2 <- sum(e^2)/length(e)
   names(sf$par)[K+1] <- "lambda"
   # names(sf$par)[K+2] <- "sigma2"
   sf$N <- dim(x)[1] 
   sf$df <- dim(x)[2] +3
   sf$loglik <- -o$val
   sf$vcov <- matrix(NA, nrow=dim(x)[2]+1, ncol=dim(x)[2]+1)
   if ( hessian == 2 || hessian == 3 ) 
      sf$vcov <- o$invhessian
   else if ( hessian == 1 )  {
      if ( !is.null(o$hessian) ) sf$vcov <- genInv(o$hessian)
   } else
   	sf$vcov <- matrix(NA, nrow=dim(x)[2]+1, ncol=dim(x)[2]+1)
   # Standard error of all parameters
   if (DEBUG & !is.null(o$hessian)) cat("Determinant for Hessian = ", det(o$hessian),"\n" )
   ## Standard error of the parameters in the production function
   ## std.err[1:3]
   # t-ratios
   if (sum(!is.na(sf$vcov))) {
	   sf$std.err <- sqrt(diag(sf$vcov))
   	sf$t.value <- sf$par/sf$std.err
   } else {
   	sf$std.err <- NA
   	sf$t.value <- NA
   }
   # t-ratio for production function
   #  (o$par/std.err)[1:(1+dim(x)[2])]

# print(loglik(o$par))
   return(sf)
} ## sfa



ebeta <- function(x,y,beta)  {
   #  Calculate residuals
   #  There is no intercept in the x matrix
   #  %*% is the inner matrix product
   y - beta[1] - x %*% beta[2:(dim(x)[2]+1)]
}



genInv <- function(X, tol = sqrt(.Machine$double.eps))
{
## Generalized Inverse of a Matrix
## Calculates the Moore-Penrose generalized inverse of a matrix X
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
    dimnames = dnx[2:1])
}



print.sfa <-  # function(x,...)  {
function (x, digits = max(3, getOption("digits") - 3), ...) 
{
#    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("\n")
    invisible(x)

#   a <- cbind("Parameters"=x$par,"  Std.err"=x$std.err,"  t-value"=x$t.value)
#   if ( length(x$t.value) > 0 ) {
#      a<-cbind(a," Pr(>|t|)"=as.integer(1000*
#	   2*pt(abs(x$t.value),x$N-x$df,lower.tail=F) )/1000)
#   }
#   print(a,digits=4,quote=F,...)
#   cat("sigma2        ",format(x$sigma2,digits=5),"\n")
#   invisible(a)
}  ## print.sfa



summary.sfa <- function(object, ...)  {
#   print.sfa(object)
   a <- cbind("Parameters"=object$par,"  Std.err"=object$std.err,
     "  t-value"=object$t.value)
   if ( length(object$t.value) > 0 ) {
      a<-cbind(a," Pr(>|t|)"=as.integer(1000*
	   2*pt(abs(object$t.value),object$N-object$df,lower.tail=F) )/1000)
   }
   print(a,digits=4,quote=F,...)
   cat("sigma2        ",format(object$sigma2,digits=5),"\n")
   cat("sigma2v = ", object$sigma2/(1 + object$lambda^2),
   ";  sigma2u = ", 
         object$sigma2*object$lambda^2/(1 + object$lambda^2),"\n")
   cat("log likelihood = ",object$loglik,"\n")
   cat("Convergence = ",object$convergence,
   	"; number of evaluations of likelihood function", object$info["neval"], "\n")
   cat("Max value of gradien:", object$info["maxgradient"], "\n")
   cat("Length of last step:", object$info["laststep"], "\n")
   cat("Final maximal allowed step length:", object$info["stepmax"], "\n")
   # print(object$count)
}  ## summary.sfa


fitted.sfa <- function(object, ...)  {
   return(object$fitted.values)
}


residuals.sfa <- function(object, ...)  {
   val <- object$residuals
   attr(val, "nobs") <- object$N
   attr(val, "df") <- object$df
   class(val) <- "residuals"
   return(val)
} ## residuals.sfa


vcov.sfa <- function(object, ...)  {
   return(object$vcov)
}


logLik.sfa <- function(object, ...)  {
   val <- -object$value
   attr(val, "nobs") <- object$N
   attr(val, "df") <- object$df
   class(val) <- "logLik"
   return(val)
} ## logLik.sfa


coef.sfa <- function(object, ...)  {
   return(object$coef)
}


eff.sfa <- function(object, type="BC", ...)  {
  switch( type,
          "BC" = return(te.sfa(object)),
          "Mode" = return(teMode.sfa(object)),
          "J" = return(teJ.sfa(object)),
          "add" = return(te.add.sfa(object)),
          "add" = return(te.add.sfa(object)),
          warning("Unknown type:", type)
        )
}
efficiencies.sfa <- eff.sfa



## Beregning af teknisk effektivitet
te.sfa <- function(object)  {
  # Hjaelpevariabler
  lambda <- object$lambda
  s2 <- object$sigma2
  ustar <- -object$residuals*lambda^2/(1+lambda^2)
  sstar <- lambda/(1+lambda^2)*sqrt(s2)

  # Teknisk efficiens for hver enhed
  TE = pnorm(ustar/sstar -sstar)/pnorm(ustar/sstar) * exp(sstar^2/2 -ustar)
  colnames(TE) <- "te"
   return(array(TE))
#   class(TE) <- "te"
#   TE
}
teBC.sfa <- te.sfa

teMode.sfa <- function(object)  {
  # Hjaelpevariabler
  lambda <- object$lambda
  ustar <- -object$residuals*lambda^2/(1+lambda^2)

  # Teknisk efficiens for hver enhed
  TE1 = matrix(exp(pmin(0,-ustar)),ncol=1)
  colnames(TE1) <- "teM"
  return(array(TE1))
}
# te1.sfa <- teMode.sfa

teJ.sfa <- function(object)  {
  # Hjaelpevariabler
  lambda <- object$lambda
  s2 <- object$sigma2
  ustar <- -object$residuals*lambda^2/(1+lambda^2)
  sstar <- lambda/(1+lambda^2)*sqrt(s2)

  # Teknisk efficiens for hver enhed
  TE2 = exp(-ustar -sstar*( dnorm(ustar/sstar)/pnorm(ustar/sstar) ) )
  colnames(TE2) <- "teJ"
  return(array(TE2))
}
# te2.sfa <- teJ.sfa


te.add.sfa <- function(object)  {
  e <- residuals.sfa(object)
  s2 <- sigma2.sfa(object)
  lambda <- lambda.sfa(object)
  # auxiliary variables
  sstar <- lambda/(1+lambda^2)*sqrt(s2)
  estar <- e * lambda / sqrt(s2)
  uJ <- sstar * (dnorm(estar)/(1 - pnorm(estar)) - estar)
  teAdd <- 1 - uJ/object$fitted.values
  class(teAdd) <- "matrix"
  colnames(teAdd) <- "teAdd"
  return(array(teAdd))
}
eff.add.sfa <- te.add.sfa


sigma2u.sfa <- function(object)  {
  s2u <- object$lambda^2 / (1+object$lambda^2) * object$sigma2
  names(s2u) <- "sigma2u"
  return(s2u)
}

sigma2v.sfa <- function(object)  {
  s2v <- object$sigma2 / (1+object$lambda^2)
  names(s2v) <- "sigma2v"
  return(s2v)
}


sigma2.sfa <- function(object) {
   s2 <- object$sigma2
  names(s2) <- "sigma2"
  return(s2) 
}


lambda.sfa <- function(object)  {
   lam <- object$lambda
   names(lam) <- "lambda"
   return(lam)
}
