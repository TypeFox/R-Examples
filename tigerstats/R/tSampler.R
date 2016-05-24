#' @title Illustrating the t-statistic

#' @description An app to explore the distribution of the t-statistic.  The user takes one sample at
#' a time from a given population.  Graphical
#' output updates the empirical distribution of the sample mean.
#' 
#' @rdname tSampler
#' @usage tSampler(form,data,max.sample.size=30,show.sample=FALSE)
#' @param form An object of class formula, of the form ~x, where x is a numeric variable from the data frame supplied by:
#' @param data A dataframe, representing the imaginary population.
#' @param max.sample.size Maximum sample size on the slider.
#' @param show.sample  If TRUE, the complete sample will be output to the console, in addition to the summary information. 
#' @return Graphical and numerical output.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @note  Uses \code{manipulate}.
#' @examples
#' \dontrun{
#' data(imagpop)
#' if (require(manipulate)) tSampler(~income,data=imagpop)
#' }
tSampler <-
function(form,data,
                        max.sample.size=30,
                        show.sample=FALSE) {
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  #We are NOT gonna monkey around with missing values:
  data <- data[complete.cases(data),] 
  dcurven <- 512
  prsd <- with(data,ParseFormula(form))
  varname <- as.character(prsd$rhs)
  numvar <- data[,varname]
  pop.mean <- mean(numvar)
  if (max.sample.size>length(numvar))  {
    stop("Sample size must be less than pop size")
  }
  exclamation <- c("Surprise!",
                   "Woo-Hoo!",
                   "Look, Ma:\nNo Hands!",
                   "I'm\nWAY\nout!",
                   "Yow!",
                   "Groovy,\nMan!",
                   "Yippee!",
                   "Far Out ...",
                   "I'm Unusual")
  t.stats <- NULL
  nprev <- NULL
  beginning <- TRUE
  dpop <- NULL
  
  manipulate(
    n=slider(2,max.sample.size,step=1,initial=2,
                label="Sample Size n"),
    SurpriseLines=checkbox(label="Show Rough Surprise Lines"),
    tcurve=checkbox(label="Show t-Curve"),
    
    {OurSample <- sample(numvar,n,replace=FALSE)
    samp.mean <- mean(OurSample)
    se.mean <- sd(numvar)/sqrt(n)
    t.stat <- (samp.mean-mean(numvar))/se.mean
     
     if (!identical(n,nprev)) {
       t.stats <<- NULL
    }
     
    if (!beginning) t.stats <<- c(t.stats,t.stat)
    
    if (show.sample  && !beginning)  {
      cat("The sample was\n")
      print(OurSample)
      cat("\n")
    }
    if (!beginning) {
    cat("t-statistic was",round(t.stat,3),"\n")
    }
    
     nprev <<- n
     if (beginning) {
       dpop <<- density(numvar,n=1024,from=min(numvar),to=max(numvar))
       plot(dpop$x,dpop$y,type="l",col="red",lwd=2,
            xlab=varname,ylab="Density",
            main="Population Distribution")
       beginning <<- FALSE
     }
     else {
       if(length(t.stats) >=2) {
       dts <- density(t.stats,n=dcurven)
       plot(dts$x,dts$y,type="l",col="blue",xlim=c(-4,4),
              xlab="x-bar",ylab="Estimated Density",
              main="Distribution of t-Statistics")
       rug(t.stats)
       points(t.stat,0,col="blue",pch=19)
       if (SurpriseLines) {
         abline(v=-2,lwd=1.5,col="red")
         abline(v=2,lwd=1.5,col="red")
         m <- length(t.stats)
         surprise.prop <- length(t.stats[abs(t.stats)>2])/m
         mtext(paste(round(100*surprise.prop,3),"% surpises so far"))
         if (abs(t.stat) > 2) {
           text(x=t.stat,y=0.1,labels=sample(exclamation,size=1))
         }
       }
       if (tcurve) {
         xvals <- seq(-4,4,by=0.01)
       points(xvals,dt(xvals,df=n-1),type="l",col="red",lwd=2)
       legend(x="topright",legend=c("t-curve","t-statistics"),col=c("red","blue"),lwd=c(2,1))
        }
       }
     }
}
  )# end manipulate
} #end tSampler

if(getRversion() >= "2.15.1")  utils::globalVariables(c("SurpriseLines","tcurve"))
