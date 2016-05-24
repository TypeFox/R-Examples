## here we assuming that the user has only one x-variabe in his fitted 
## model and that this x-variable may have been transformed into a power
## x^power before fitting
## so for predictive values we need the new x values (xvalues) 
## the original name of the x variable (xname)
## and possible its power (power).
## Mikis Stasinopoulos Monday, August 2, 2004
centiles.pred <- function(obj,
                        type = c("centiles", "z-scores", "standard-centiles"), 
                       xname = NULL,
                     xvalues = NULL,   
                       power = NULL,
                        yval = NULL,
                        cent = c(.4,2,10,25,50,75,90,98,99.6),
                         dev = c(-4,-3,-2,-1,0,1,2,3,4),
                        plot = FALSE,
                      legend = TRUE, 
                              ...)
{
## ------function calc.cent-----------------------------------------------------
# Huiqi  
 calc.cent <- function(xvar, cent) {
         o <- order(xvar)
       mat <- xvar[o]
      cent <- cent
        for (var in cent) {
            if (lpar == 1) {        
                newcall <- call(qfun, var/100, mu = mu[o])
            }
            else if (lpar == 2) {
                newcall <- call(qfun, var/100, mu = mu[o], 
                  sigma = sigma[o])
            }
            else if (lpar == 3) {
                newcall <- call(qfun, var/100, mu = mu[o], 
                  sigma = sigma[o], nu = nu[o])
            }
            else {
                newcall <- call(qfun, var/100, mu = mu[o], 
                  sigma = sigma[o], nu = nu[o], 
                  tau = tau[o])
            }
            ll <- eval(newcall)
            mat <- cbind(mat, ll)
        }
        mat <- as.data.frame(mat)
        nnn <- paste("C", as.character(cent), sep = "")
        names(mat) <- c(xname, nnn)
        return(mat)
    }
 # calc.cent <- function(xvar, cent)
 #  {
 #       mat <- xvar
 #      cent <- cent          
 #     for(var in cent) 
 #       { 
 #         if(lpar==1) 
 #         {
 #         newcall <-call(qfun,var/100,
 #                   mu=mu[order(xvar)]) 
 #         }
 #        else if(lpar==2)
 #         {
 #         newcall <-call(qfun,var/100,
 #                   mu=mu[order(xvar)],
 #                   sigma=sigma[order(xvar)]) 
 #         }
 #        else if(lpar==3)
 #         { 
 #         newcall <-call(qfun,var/100,
 #                   mu=mu[order(xvar)],
 #                   sigma=sigma[order(xvar)],
 #                   nu=nu[order(xvar)])
 #         }
 #       else 
 #         {
 #         newcall <-call(qfun,var/100,
 #                   mu=mu[order(xvar)],
 #                   sigma=sigma[order(xvar)],
 #                   nu=nu[order(xvar)],
 #                   tau=tau[order(xvar)]) 
 #         }
 #         ll <- eval(newcall)
 #        mat <- cbind(mat,ll)
 #      }
 #         mat <- as.data.frame(mat)
 #         nnn <- paste("C", as.character(cent), sep="")
 #  names(mat) <-c(xname,nnn)
 #  return(mat)
 #  }  
## end of function calc.cent
   plot.mat <- function(mat, cent, legend ,...)
   {
    lcent <- dim(mat)[2]
     xleg <- min(mat[,1])
     yleg <- max(mat[,2:lcent])
    plot(mat[,1],mat[,2],ylim=c(min(mat[,2:lcent]),max(mat[,2:lcent])), type="n",...)
     for(i in 2:lcent)    lines(mat[,1],mat[,i],col=i)
     if (legend)
      legend(list(x=xleg,y=yleg), legend = cent, col=c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), lty=1, ncol=1, bg="white")# 
     invisible()
   }
##-----------------------------------------------------
## the main function start here
##  checking the object
if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
## checking the xvalues
if (is.null(xvalues)) stop(paste("The xvalues  argument is not specified", "\n", ""))
## checking the xname
if (is.null(xname)) stop(paste("The xname argument is not specified", "\n", ""))
if (!is.character(xname)) stop(paste("The xname argument is not a character", "\n", ""))
## checking for continuous family
#if(!obj$type=="Continuous") 
#   stop(paste("The centiles are working only with continuous distributions", "\n", ""))
## if power 
    xvar <- if (!is.null(power))   xvar <-  xvalues^power
           else xvalues
## create a data frame
    newx <- data.frame(xvar)
colnames(newx) <- xname 
    lpar <- length(obj$parameters)
## the problem here is that if any parameters is fixed the predict will not work
## here is the fix MS Wednesday, June 28, 2006 
if ("mu"%in%obj$parameters )
    {if ( is.null(obj$mu.fix))    
      mu <- predict(obj,what = "mu",   newdata = newx, type = "response", ...)
  else if (obj$mu.fix==TRUE) mu <- rep(fitted(obj, "mu")[1], length(xvar)) 
  }
if ("sigma"%in%obj$parameters)
   {if (is.null(obj$sigma.fix))    
      sigma <- predict(obj,what = "sigma",   newdata = newx, type = "response", ...)
  else if (obj$sigma.fix==TRUE) sigma <- rep(fitted(obj, "sigma")[1], length(xvar)) 
  }
if ("nu"%in%obj$parameters )
   { if  (is.null(obj$nu.fix))    
      nu <- predict(obj,what = "nu",   newdata = newx, type = "response", ...)
  else if (obj$nu.fix==TRUE) nu <- rep(fitted(obj, "nu")[1], length(xvar))
   } 
if ("tau"%in%obj$parameters )
   { if (is.null(obj$tau.fix))    
      tau <- predict(obj,what = "tau",   newdata = newx, type = "response", ...)
  else if (obj$tau.fix==TRUE) tau <- rep(fitted(obj, "tau")[1], length(xvar)) 
  }  
## which type
    type <- match.arg(type)
## only for type centiles
if (type=="centiles")
     {   
    fname <- obj$family[1]
     qfun <- paste("q",fname,sep="")
     #  Title <- paste("Centile curves using",fname, sep=" ")
     xvar <- xvalues
   # oxvar <- xvar[order(xvar)]
      mat <- calc.cent(xvar=xvar, cent=cent)
      if (plot)
      plot.mat(mat,cent,legend,...)
      return(mat)
     }
##   only for type z-scores  
if (type=="z-scores")
     {
      if (is.null(yval)) stop("the y values should be set if type=z-scores is used")
      if (length(yval)!= length(xvalues)) stop("length of xvalues and yval is not the same")
    fname <- obj$family[1]
     qfun <- paste("p",fname,sep="")
         if(lpar==1) 
          {newcall <-call(qfun,yval,mu=mu) }
         else if(lpar==2)
          {newcall <-call(qfun,yval,mu=mu,sigma=sigma) }
         else if(lpar==3)
          {newcall <-call(qfun,yval,mu=mu,sigma=sigma,nu=nu) }
        else 
          {newcall <-call(qfun,yval,mu=mu,sigma=sigma,nu=nu,tau=tau) }
          cdf <- eval(newcall)       
      rqres <- qnorm(cdf)
      return(rqres) 
     }
## only for type standard-centiles
if (type=="standard-centiles")
     {   
     cent<-pnorm(dev)*100
    fname <- obj$family[1]
     qfun <- paste("q",fname,sep="")
     #  Title <- paste("Centile curves using",fname, sep=" ")
     xvar <- xvalues
      mat <- calc.cent(xvar=xvar, cent=cent)
      nnn <- paste(as.character(dev), sep="")
   names(mat) <-c(xname,nnn)
    if (plot)
      plot.mat(mat,dev,legend,...)
      return(mat)
     }        

} 
