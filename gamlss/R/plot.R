# this is the plot.gamlss function 
# created by PA  May 2002
# last change by MS Friday, Wednesday, December 17, 2003 at 09:11
# to incoorporate options in the plotying parameters
# the following options have been used for the BCT paper 
# par(mfrow=c(2,2), mar=par("mar")+c(0,1,0,0), col.axis="blue4", col="blue4", col.main="blue4",col.lab="blue4",pch="+",cex=.45, cex.lab=1.2, cex.axis=1, cex.main=1.2  )
# it fails in  

#plot(mod)
#Error in xy.coords(x, y) :
#'x' is a list, but does not have components 'x' and 'y'
#Calls: plot ... axis -> points -> points.default -> plot.xy -> xy.coords
#Execution halted
# plot ... axis -> points -> points.default -> plot.xy -> xy.coords
# brian
# Error in xy.coords(x, y) :
#  'x' is a list, but does not have components 'x' and 'y'
#Calls: plot ... axis -> points -> points.default -> plot.xy -> xy.coords#

#This is less obvious: the call is

#         rug(residx, col="red", points(par(col="blue4")))

#and the baffling call to points() is the problem.

#-------------------------------------------------------------------------------
plot.gamlss <- function (x, xvar=NULL, parameters=NULL, ts=FALSE, summaries=TRUE, ...) 
{
    if (!is.gamlss(x))  stop(paste("This is not an gamlss object", "\n", ""))
## chech for the residuals 
    if (is.null(x$residuals)) #
         stop(paste("There are no quantile residuals in the object"))
    residx <- resid(x) # get the residuals 
         w <- x$weights
    xlabel <- if(!missing(xvar)) deparse(substitute(xvar)) else deparse(substitute(index))
## plotting parameters
    if(is.null(parameters))
          op <- par(mfrow=c(2,2), mar=par("mar")+c(0,1,0,0), col.axis="blue4", col.main="blue4", col.lab="blue4",  col="darkgreen", bg="beige" )
    else  op <- parameters
## now the two top  figures 
## if time series plot acf and pacf  
    if(identical(ts, TRUE))
     {  # get the acf and pacf
     acf.new<-acf(residx,plot=FALSE)
     plot(acf.new,xlim=c(2,length(acf.new$acf)),ylim=range(acf.new$acf[-1]))   # ms Tuesday, August 19, 2003 at 11:04
     pacf(residx)
     }
     else 
     {# otherwise 
     ## I am assuming that is x$noObs!=x$N then we have weights (with frequencies)
     if (length(residx)==x$N)
        {
         fittedvalues <- if(is.null(fitted(x))) fitted(x,"sigma") else fitted(x) # MS Wednesday, September 10, 2003 at 21:20
         ## whether index or x-variable
         if(is.null(xvar))     xvar <- seq(1,length(residx),1) # MS
        }
        else
        { # if weights
         fittedvalues <- rep( if(is.null(fitted(x))) fitted(x,"sigma") else fitted(x), w)
          xvar <- if(is.null(xvar))  seq(1,length(residx),1) else rep(xvar,w)
        } 
    # top left
    plot(fittedvalues , residx,
         xlab = "Fitted Values",  
         ylab = "Quantile Residuals", 
         main = "Against Fitted Values",
         frame.plot = TRUE) 
    # top right  
    plot(xvar, residx, 
         ylab = "Quantile Residuals",
         xlab = xlabel, 
         main = paste("Against ", xlabel), 
         frame.plot = TRUE) #  points(par(col="blue4"))
     }    
    plot(density(residx), 
         xlab = "Quantile. Residuals", 
         ylab = "Density", 
         main = "Density Estimate",
         frame.plot = TRUE, 
         col="black", 
         lwd=0.4 ) #col="deepskyblue4", col="darkgreen", 
         rug(residx, col="red")
 
    qqnorm(residx, main = "Normal Q-Q Plot",
            xlab = "Theoretical Quantiles",
            ylab = "Sample Quantiles", 
            plot.it = TRUE, 
            frame.plot = TRUE, 
            col="darkgreen")
     lines(residx, residx, col="red" , lwd=.4, cex=.4 )
 
     if ( identical(summaries, TRUE))
               { 
                     qq <- as.data.frame(qqnorm(residx, plot = FALSE))
               Filliben <- cor(qq$y,qq$x)
                    # mr <- as.matrix(residx)
                    m.1 <- mean(residx)
                    m.2 <- var(residx) # cov.wt(mr,w)$cov
                  n.obs <- sum(w) 
                    m.3 <- sum((residx-m.1)**3)/n.obs 
                    m.4 <- sum((residx-m.1)**4)/n.obs 
                    b.1 <- m.3^2/m.2^3
                sqrtb.1 <- sign(m.3)*sqrt(abs(b.1))
                    b.2 <- m.4/m.2^2 
                     cat("*******************************************************************")
                     cat("\n")
                     if (identical(x$type,"Continuous")) 
                         {cat("\t","     Summary of the Quantile Residuals")}
                     else{cat("\t","Summary of the Randomised Quantile Residuals")}    
                     cat("\n")
                     cat("                           mean   = ", m.1, "\n")
                     cat("                       variance   = ", m.2, "\n")
                     cat("               coef. of skewness  = ", sqrtb.1, "\n")
                     cat("               coef. of kurtosis  = ", b.2, "\n")
                     cat("Filliben correlation coefficient  = ", Filliben, "\n")
                     cat("*******************************************************************")
                     cat("\n")

               }
    par(op)
}
#par(mfrow=c(2,2), mar=par("mar")+c(0,1,0,0), col.axis="blue4", col="blue4", col.main="blue4",col.lab="blue4",pch="+",cex=.45, cex.lab=1.2, cex.axis=1, cex.main=1.2  )
#---------------------------------------------------------------------------------------
