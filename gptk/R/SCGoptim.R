SCGoptim <-
function (x, fn, grad, options, ...) {
  ## options = list(maxit, ln, xtol, fnTol, optimiser="SCG", gradcheck=FALSE)
  ##cat ("\n SCG Optimisation begins! \n")

  if ( "maxit" %in% names(options) && !is.na(options$maxit) ) {
    niters <- options$maxit
  } else {
    niters <- 100
  }

  display <- options$display
  gradcheck <- options$gradcheck
  
  ## y = fn (x)
  func <- function(x, ...) fn(x, ...)
      
  ## gradient function = gr (x)
  gradfunc <- function(x, ...) grad(x, ...)

  nparams <- length(x)

  sigma0 <- 1e-4
  fold <- func(x, ...)
  fnow <- fold
  gradnew <- gradfunc(x, ...)
  gradold <- gradnew
  d <- -gradnew
  success <- 1
  nsuccess <- 0
  beta <- 1
  betamin <- 1e-15
  betamax <- 1e100
  eps <- 2.2204e-16
  j <- 1

  while ( j<=niters ) {
    if ( success == 1 ) {
      mu <- crossprod(d, gradnew)
      if ( mu>=0 ) {
        d <- -gradnew
        mu <- crossprod(d, gradnew)
      }
      kappa <- crossprod(d, d)
      if ( kappa<eps ) {
        xmin <- x
        objective <- fnow
        ans <- list(xmin=xmin, objective=objective)
        return (ans)
      }
      sigma <- (sigma0/sqrt(kappa))[1]
      xplus <- x+sigma*d
      gplus <- gradfunc(xplus, ...)
      theta <- crossprod(d, gplus-gradnew)/sigma
    }

    delta <- theta + beta*kappa
    if ( delta<=0 ) {
      delta <- beta*kappa
      beta <- beta - theta/kappa
    }
    alpha <- (-mu/delta)[1]

    xnew <- x+alpha*d
    fnew <- try( func(xnew, ...), silent=TRUE )
    if ( !is.finite(fnew) ) fi <- 1
    while ( !is.finite(fnew) ) {
      if (display) {
        if ( fi==1 ) {   
	  message("\t function evaluation failed in SCG.")      
        } else {      
          message(".")
        }
      }
      alpha <- alpha/2
      xnew <- x+alpha*d
      fnew <- try( func(xnew, ...) )
      fi <- fi+1
      if ( is.finite(fnew) ) {
        if (display)
          message("\n")
        fi <- 0
      }
    }

    Delta <- 2*(fnew-fold)/(alpha*mu)
    if ( Delta>=0 ) {
      success <- 1
      nsuccess <- nsuccess+1
      x <- xnew
      fnow <- fnew
    } else {
      success <- 0
      fnow <- fold
    }
    if (display)
      cat("Cycle ", j, "Error ", round(fnow, digits=4), "Scale ", beta, "\n")

    if ( success == 1 )
      if ( (max(abs(alpha*d))<options$xtol) & (max(abs(fnew-fold))<options$fnTol) ) {
        xmin <- x
        objective <- fnew
        ans <- list(xmin=xmin, objective=objective)
        return (ans)
      } else {
        fold <- fnew
        gradold <- gradnew
        gradnew <- gradfunc(x, ...)
        if ( crossprod(gradnew, gradnew) == 0 ) {
          xmin <- x
          objective <- fnew
          ans <- list(xmin=xmin, objective=objective)
          return (ans)
        }
      }

    if ( Delta < 0.25 ) 
      beta <- min(2*beta, betamax)

    if ( Delta > 0.75 )
      beta <- max(0.5*beta, betamin)

    if ( nsuccess == nparams ) {
      d <- -gradnew
      nsuccess <- 0
    } else {
      if ( success==1 ) {
        gamma <- (crossprod(gradold-gradnew, gradnew)/mu)[1]
        d <- gamma*d-gradnew
      }
    }

    j <- j+1
  }

  xmin <- x
  objective <- fold
  ans <- list(xmin=xmin, objective=objective)
  warning("Maximum number of iterations has been exceeded.\n")
  return (ans)
}
