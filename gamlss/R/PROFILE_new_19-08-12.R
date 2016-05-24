prof.dev<- function ( object, 
                      which = NULL,  
                        min = NULL, 
                        max = NULL, 
                       step = NULL,
                     length = 7, 
               startlastfit = TRUE, 
                       plot = TRUE,
                       perc = 95, 
                      ...) 
{
    cat("*******************************************************************", "\n")
   if (is.null(object$y) ) 
       stop(paste(deparse(substitute(object)), "does not have `y' components"))
   if (is.null(min)) 
       stop(paste("The min has not been specified"))
   if (is.null(max)) 
       stop(paste("The max has not been specified"))
 #  if (is.null(step)) 
 #      stop(paste("The step has not been specified"))
   if (is.null(which)) 
       stop(paste("The parameter has not been specified"))
   else 
    {
       if (!any(which == c("mu", "sigma", "nu", "tau")))  
           stop(paste("The parameter has not been specified correctly."))
    }

    #interval <- seq(from=min, to=max, by=step)
         interval <- if (is.null(step)) seq(from=min, to=max, length.out=length) else seq(from=min, to=max, by=step)
         interval <- ifelse(interval==0,0.0001,interval)
    fix.parameter <- rep(0, length(interval))
      G.deviances <- rep(0, length(interval))
                m <- length(interval)
    for (i in 1:m)
    {
        if (which == "mu") 
        {   
            if ("mu"%in%object$parameters) 
            {
                        mu.start <- rep(interval[i],length(object$y))
              object$call$mu.fix <- TRUE
            object$call$mu.start <- mu.start
                cat("mu.start=(" ,interval[i],")","\n")
            }
            else
                stop(paste("The mu is not a valid parameter of the current model."))
        }
        if (which == "sigma")
        {
            if ("sigma"%in%object$parameters)  
            {
                      sigma.start <- rep(interval[i],length(object$y))
            object$call$sigma.fix <- TRUE
          object$call$sigma.start <- sigma.start
                cat("sigma.start=(" ,interval[i],")","\n")
            }
            else
                stop(paste("The sigma is not a valid parameter of the current model."))
        }
        if (which == "nu")
        {
            if ("nu"%in%object$parameters) 
            {
                         nu.start <- rep(interval[i],length(object$y))
               object$call$nu.fix <- TRUE
             object$call$nu.start <- nu.start
                cat("nu.start=(" ,interval[i],")","\n")
            }
            else
                stop(paste("The nu is not a valid parameter of the current model."))
        }
        if (which == "tau") 
        {
            if ("tau"%in%object$parameters)
            {
                         tau.start <- rep(interval[i],length(object$y))
               object$call$tau.fix <- TRUE
             object$call$tau.start <- tau.start
                cat("tau.start=(" ,interval[i],")","\n")
            }
            else
                stop(paste("The tau is not a valid parameter of the current model."))
        }
        
        #                    FIT THE MODEL     
        object$call$control <- object$control
                        out <- eval(object$call)
           fix.parameter[i] <- interval[i]
             G.deviances[i] <- out$G.deviance
            if (startlastfit)
        {
            if ("mu"%in%object$parameters)       object$call$mu.start <- out$mu.fv
            if ("sigma"%in%object$parameters) object$call$sigma.start <- out$sigma.fv
            if ("nu"%in%object$parameters)       object$call$nu.start <- out$nu.fv
            if ("tau"%in%object$parameters)     object$call$tau.start <- out$tau.fv
        }
        cat("*******************************************************************", "\n")
    }
 
 prof.out <- cbind(interval, G.deviances)
       mx <- which.min(G.deviances)
       if ((mx<=1)||(mx>=m)) stop("increase the interval to contain the MLE of the parameter") 
 prof.fun <- splinefun(interval, G.deviances) # the profile likelihood as a function
     PML <- uniroot(prof.fun, c(min, max), deriv=1)$root # get the ML
    Gmin <- prof.fun(PML)      
     lim <- Gmin + qchisq((perc/100), 1)
      xl <- as.vector(interval) 
      CI <- c(NA, NA)
if (plot)
  { 
   #
    loglik <- G.deviances
    xlabel <- paste("Grid of the",which,"parameter")
      xlab <- xlabel
      ylab <- "Global Deviances"
      main <- "Profile Global Deviance"
     curve(prof.fun, min, max, xlab=xlab, ylab=ylab, main=main, col="darkgreen", frame.plot = TRUE)
     plims <- par("usr")  
     segments(PML, plims[3], PML, Gmin, lty = 3)
        la <- xl[mx]
        y0 <- plims[3]
    if (lim < max(loglik))
    {
        abline(h = lim, lty = 3)
        y0 <- plims[3]
        scal <- (1/10 * (plims[4] - y0))/par("pin")[2] #par("pin")[2]=the height of the plot in inches
        scx <- (2/10 * (plims[2] - plims[1]))/par("pin")[1] #par("pin")[1]=the width of the plot in inches
        # MS change to 2/10 from 1/10 Sunday, December 9, 2007 at 23:29
        text(xl[1] + scx, lim + scal, paste(perc,"%"))
     }
#    #Defines the lower bound
    if (loglik[1] > lim) {
                 leftFun <- function(x)  {prof.fun(x)-lim} # 
               lcrossing <- uniroot(leftFun, c(min, PML))$root # get the ML
                   CI[1] <- lcrossing
                   segments( lcrossing, y0,  lcrossing, lim, lty = 3)
    }
    #Defines the upper bound
    if (loglik[m] > lim) {
               rcrossing <- uniroot(leftFun, c(PML, max))$root # get the ML
              segments( rcrossing, y0,  rcrossing, lim, lty = 3)
               CI[2] <- rcrossing
    }
    cat("*******************************************************************", "\n")
    cat("The Maximum Likelihood estimator is " ,PML, "\n")
    cat("with a Global Deviance equal to ", Gmin, "\n")
    if ((loglik[1] > lim) && (loglik[m] > lim))    
        {cat("A ", perc,"% Confidence interval is: (" ,lcrossing, ",", rcrossing, ") \n")}
    cat("*******************************************************************", "\n")           
          out <-list(values=prof.out, fun=prof.fun, min=min, max=max, max.value=PML, CI=CI, criterion="GD" )
        class(out) <- "ProfLikelihood.gamlss"
     return(invisible(out))
  }
    else return(prof.out)
}
#prof.dev1(h,"nu",min=-2.000,length=6,max=2,type="l")