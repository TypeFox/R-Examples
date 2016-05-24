`print.MCepi` <-
function(x, ...)
  {
    s <- summary(x, ...)
    print(s, ...)
  }

summary.MCepi <- function(object, ...)
  {
    rl <- list(obj=object)    
    class(rl) <- "summary.MCepi"
    
    rl$final <- data.frame(rbind(getvac(object), getcost(object)),
                           row.names=c("vac", "cost"))
    
    rl
  }

print.summary.MCepi <- function(x, ...)
  {
    ## print information about the call
    cat("\nCall:\n")
    print(x$obj$call)

    ## print information on the final number(s)
    ## of vaccinated people
    cat("\nDistribution of vaccinations administered:\n")
    print(x$final[1,])

    ## print information on the final cost
    cat("\nDistribution of final costs:\n")
    print(x$final[2,])

    cat("\n")

  }

summary.epiman <- function(object, ...)
  {
    rl <- list(obj=object)    
    class(rl) <- "summary.epiman"

    rl$params <- summary(object$samp)

    rl$final <- data.frame(vac=getvac(object), cost=getcost(object))

    rl
  }

print.summary.epiman <- function(x, ...)
  {
    ## print information about the call
    cat("\nCall:\n")
    print(x$obj$call)

    ## print information on the estimates of the parameters
    ## in the final time step
    cat("\nDistribution of SIR model parameters:\n")
    print(x$params)

    ## print infomrmation on the final cost and vacs
    cat("\nVaccinations administered and cost:\n")
    print(x$final)

    cat("\n")
  }

print.epiman <- function(x, ...)
  {
    s <- summary(x, ...)
    print(s, ...)
  }


summary.optvac <- function(object, ...)
  {
    rl <- list(obj=object)    
    class(rl) <- "summary.optvac"

    rl$best <- getpolicy(object)
    rl$worst <- getpolicy(object, "worst")

    rl
  }

print.summary.optvac <- function(x, ...)
  {
    ## print information about the call
    cat("\nCall:\n")
    print(x$obj$call)

    cat("\nThe cost matrix is contained in the $C field.\n")

    cat("\nThe best and worst policies are:\n")
    print(rbind(x$best, x$worst))

    cat("\n")
  }

print.optvac <- function(x, ...)
  {
    s <- summary(x)
    print(s, ...)
  }
