#' @title Distribution of the Sample Mean

#' @description An app to explore the Central Limit Theorem.
#' 
#' @rdname SampDistMean
#' @usage SampDistMean(pop,max.samp.size=50,sim.reps=1000)
#' @param pop A data frame representing the population from which samples are taken.
#' @param max.samp.size  Largest sample size shown on the slider.
#' @param sim.reps Number of simulation repetitions to construct empirical distribution of the sample mean.
#' @return Graphical and numerical output.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @note  Uses \code{manipulate} in RStudio.
#' @examples
#' \dontrun{
#' data(imagpop)
#' if (require(manipulate)) SampDistMean(imagpop)
#' }
SampDistMean <- function(pop,max.samp.size=50,sim.reps=1000) {
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  #pop should be a data frame
  
  #We are NOT gonna monkey around with missing values:
  pop <- pop[complete.cases(pop),]
  
  #remove non-numerical variables from pop
  numvars <- NULL
  for(name in names(pop))  {
    numvars <- c(numvars,is.numeric(pop[,name]))
  }
  pop <- pop[,which(numvars)]
  
  manipulate(
    n=slider(1,max.samp.size,initial=1,label="Sample Size n"),
    variable=picker(as.list(names(pop))),
    curvetype=picker("None","Density Estimate","Theoretical Normal"),
{varname <- as.character(variable)
 popvals <- pop[,varname]
 popmin <- min(popvals)
 popmax <- max(popvals)
 
 #determine upper y limit for histograms
 ymax <- 1.2*sqrt(max.samp.size)/(sqrt(2*3.1416)*sd(popvals))
 
 hist(popvals,freq=FALSE,ylim=c(0,ymax),
      main=paste("Histogram of",varname),
      xlab=varname,col="red")
 
 reps <- sim.reps
 sim.means <- replicate(reps,mean(sample(popvals,n,replace=FALSE)))
 hist(sim.means,col=rgb(0, 0, 1,0.5),freq=FALSE,add=TRUE)
 legend("topright", c("Population",paste(reps,"Sample Means")),
        fill = c("red", rgb(0, 0, 1,0.5)), title="Distributions",
        cex=0.75)
 
 if(curvetype=="Density Estimate") {
   d <- density(sim.means)
   lines(d$x,d$y,col="red")}
 
 if (curvetype=="Theoretical Normal")  {
   popsize <- length(popvals)
   sdmean <- sd(popvals)/sqrt(n)*sqrt((popsize-n)/(popsize-1))
   popmean <- mean(popvals)
   curve(dnorm(x,mean=popmean,sd=sdmean),add=TRUE,col="red",n=1001)
 }
 
 
 quantities <- c("population mean mu","expected value of sample mean",
                 "st dev of pop","st dev of sample mean")
 theoretical <-c(mean(popvals),mean(popvals),sd(popvals),sd(popvals/sqrt(n)))
 in.simulation <- c(NA,mean(sim.means),NA,sd(sim.means))
 frm <- data.frame(theoretical,in.simulation)
 rownames(frm) <- quantities
 print(frm)
}
  )
}  #end of SampDistMeans

if(getRversion() >= "2.15.1")  utils::globalVariables(c("picker", "variable","curvetype"))
