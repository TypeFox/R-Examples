#' @title Repeated Sampling for a Mean (Slow)

#' @description An app to explore the sampling distribution of the sample mean.  The user takes one sample at
#' a time from a given population.  Output to the console describes relevant features of the sample, and graphical
#' output updates the empirical distribution of the sample mean.
#' 
#' @rdname MeanSampler
#' @usage MeanSampler(form,data,max.sample.size=30,show.sample=FALSE)
#' @param form an object of class formula, of the form ~x, where x is a numeric variable from the data frame supplied by:
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
#' if (require(manipulate)) MeanSampler(~income,data=imagpop)
#' }
MeanSampler <-
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
  
  sampmeans <- NULL
  nprev <- NULL
  beginning <- TRUE
  dpop <- NULL
  
  manipulate(
    n=slider(1,max.sample.size,step=1,initial=1,
                label="Sample Size n"),
    {OurSample <- sample(numvar,n,replace=FALSE)
    samp.mean <- mean(OurSample)
     
     if (!identical(n,nprev)) {
       sampmeans <<- NULL
    }
     
    if (!beginning) sampmeans <<- c(sampmeans,samp.mean)
    
    if (show.sample  && !beginning)  {
      cat("The sample was\n")
      print(OurSample)
    }
    if (!beginning) {
    cat("x-bar was",round(samp.mean,3),"\n")
    cat("mu is",round(pop.mean,3),"\n")
    }
    
     nprev <<- n
     if (beginning) {
       dpop <<- density(numvar,n=1024,from=min(numvar),to=max(numvar))
       plot(dpop$x,dpop$y,type="l",col="red",lwd=2,
            xlab=varname,ylab="Density",
            main="Population Distribution")
     }
     else {
       if(length(sampmeans) >=2) {
       dsamp <- density(sampmeans,n=dcurven,from=min(numvar),to=max(numvar))
       plot(dsamp$x,dsamp$y,type="l",col="blue",
              xlab="x-bar",ylab="Estimated Density",
              main="Distribution of Sample Means")
       points(dpop$x,dpop$y,type="l",col="red",lwd=2)
       rug(sampmeans)
       points(samp.mean,0,col="blue",pch=19)
       abline(v=pop.mean,lwd=1.5,col="red")
       legend(x="topright",legend=c("Population","Sample Means"),col=c("red","blue"),lwd=c(2,1))
        }
     }
     
     beginning <<- FALSE
    }
  )
}
