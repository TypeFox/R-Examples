# this is a general profile (deviance or GAIC) function 
# it can be used for profiling degres of freedom in a GAM model  
# or in fact any coefficient of a term in a GAMLSS models
# created by MS Thursday, June 19, 2003 at 13:29 
# modified by MS Monday, August 25, 2003 at 13:25 to include intervals for profile GD
# tested again and modified by MS at 16-08-12
# new features are that now there is no need to fit a lot of points since the 
# profile deviance (or GAIC) function is approximated by a spline approximation 
# also the reported maximum is more likely to be correct since ir is calculated
# from the approximate function rather that a grip of values which can be crude
prof.term <- function (model = NULL, 
                   criterion = c("GD","GAIC"), 
                     penalty = 2.5, 
                       other = NULL,
                         min = NULL, 
                         max = NULL, 
                        step = NULL, 
                      length = 7,
                      xlabel = NULL,
                        plot = TRUE,
                        perc = 95,
                  start.prev = TRUE, 
                        ...) 
{
if(is.null(model)) stop("you have not defined the model")
if(is.null(min)) stop("you have not defined the minimum")
if(is.null(max)) stop("you have not defined the maximum")
#if(is.null(step)) stop("you have not defined the step")
criterion  <- match.arg(criterion)
#if(!criterion%in%c("IC","GD")) stop("criterion should be IC or GD")
  interval <- if (is.null(step)) seq(from=min, to=max, length.out=length) else seq(from=min, to=max, by=step)
     I.C <- rep(0, length(interval)) 
    call <- model
# if (!is.null(model$data)) 
#        {
#        DaTaa<-model$data  
#        attach(eval(substitute(DaTaa)))
#        on.exit(detach(eval(substitute(DaTaa))))
#        }
    for (i in 1:length(interval))
    {
    this<<- this <- interval[i] # mikis Thursday, March 27, 2008 
    if (!is.null(other)) eval(other)
      mod.1 <- eval(call) 
       call <- mod.1$call
    if (start.prev)
    {
      if (   "mu"%in%mod.1$parameters)  call$mu.start<-fitted(mod.1,"mu")
      if ("sigma"%in%mod.1$parameters)  call$sigma.start <- fitted(mod.1,"sigma")
      if (   "nu"%in%mod.1$parameters)  call$nu.start    <- fitted(mod.1,"nu")
      if (  "tau"%in%mod.1$parameters)  call$tau.start   <- fitted(mod.1,"tau")
    }
 I.C[i]<-  if(criterion=="GD")  deviance(mod.1) else  GAIC(mod.1,k=penalty)
    } # finish the loop
      xlab <- if(!is.null(xlabel)) xlabel else "parameter" 
      ylab <- if(criterion=="GD") "Global Deviances" else   paste("GAIC pen=",penalty)
      main <- if(criterion=="GD") "Profile Global Deviance" else  "Profile GAIC"
  prof.out <- cbind(interval, I.C)
        mx <- which.min(I.C)
         m <- length(interval)
  if ((mx<=1)||(mx>=m)) stop("increase the interval to contain the MLE of the parameter")
  prof.fun <- splinefun(interval, I.C) # the profile likelihood as a function
       PML <- uniroot(prof.fun, c(min, max), deriv=1)$root # get the ML
      Gmin <- prof.fun(PML) # get the deviance
       lim <- Gmin + qchisq((perc/100), 1)
       xl <- as.vector(interval)
       CI <- c(NA, NA)
if (plot) 
   {
  # start plotting
  # first the curve
  curve(prof.fun, min, max, xlab=xlab, ylab=ylab, main=main, col="darkgreen", frame.plot = TRUE)
  plims <- par("usr")
  # then the minimum
  segments(PML, plims[3], PML, Gmin, lty = 3)
  # we need more inthe plot if GD
    if(criterion=="GD")
      { 
            if (lim < max(I.C)) 
               {
                 abline(h = lim, lty = 3)
                  y0 <- plims[3]
                  scal <- (1/10 * (plims[4] - y0))/par("pin")[2] #par("pin")[2]=the height of the plot in inches
                  scx <- (2/10 * (plims[2] - plims[1]))/par("pin")[1] #par("pin")[1]=the width of the plot in inches 
                  # MS change to 2/10, Sunday, December 9, 2007 at 23:28     
                  text(xl[1] + scx, lim + scal, paste(perc,"%") )        
                }
            #    ind <- range((1:m)[I.C < lim]) #gets the x-values for the given range that their GD is in the range GD-or+3.84
             if (I.C[1] > lim)  #Defines the lower bound 
               {
                 leftFun <- function(x)  {prof.fun(x)-lim} # 
               lcrossing <- uniroot(leftFun, c(min, PML))$root # get the ML
                   CI[1] <- lcrossing
                  segments( lcrossing, y0,  lcrossing, lim, lty = 3)
                }
            if (I.C[m] > lim)  #Defines the upper bound
               {
              rcrossing <- uniroot(leftFun, c(PML, max))$root # get the ML
              segments( rcrossing, y0,  rcrossing, lim, lty = 3)
                  CI[2] <- rcrossing
                }
                cat("*******************************************************************", "\n")
                cat("The Maximum Likelihood estimator is " ,PML, "\n")
                cat("with a Global Deviance equal to ", Gmin, "\n")
                if ((I.C[1] > lim) && (I.C[m] > lim))    
                {cat("A ", perc,"% Confidence interval is: (" ,lcrossing, ",", rcrossing, ") \n")}
                cat("*******************************************************************", "\n")             
      } # end if GD 
    else
    {
      cat("*******************************************************************", "\n")
      cat("The Mimimum is " ,PML, "\n")
      cat("with an an GAIC(",penalty,") =", Gmin,  "\n")
      cat("*******************************************************************", "\n")
    }
    } # end if plot 
        out <-list(values=prof.out, fun=prof.fun, min=min, max=max, max.value=PML, CI=if(criterion=="GD") CI else NULL, criterion=criterion)
        class(out) <- "ProfLikelihood.gamlss"
     return(invisible(out))
}