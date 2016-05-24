#' @title Distribution of the Difference of Sample Means

#' @description An app to explore the Central Limit Theorem in the context of the difference of sample means.
#' 
#' @rdname SampDist2Means
#' @usage SampDist2Means(pop,max.samp.sizes=50,sim.reps=1000)
#' @param pop A data frame representing the population from which samples are taken.
#' @param max.samp.sizes  Largest sample sizes shown on the sliders.
#' @param sim.reps Number of simulation repetitions to construct empirical distribution of difference of sample
#' means.
#' @return Graphical and numerical output.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @note  Uses \code{manipulate} in RStudio.  Also requires package \code{lattice}.
#' @examples
#' \dontrun{
#' data(imagpop)
#' if (require(manipulate)) SampDist2Means(imagpop)
#' }
SampDist2Means <- function(pop,max.samp.sizes=50,sim.reps=1000) {
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  #pop should be a data frame with numerical and factor variables
  
  #We are NOT gonna monkey around with missing values:
  pop <- pop[complete.cases(pop),]
  
  #Separate numerical variables and factors, keep only the factors having two levels
  numvars <- NULL
  for(name in names(pop))  {
    numvars <- c(numvars,is.numeric(pop[,name]))
  }
  npop <- pop[,which(numvars)]
  
  fpop <- pop[,which(!numvars)]
  f2vars <- NULL
  for (name in names(fpop))  {
    f2vars <- c(f2vars,nlevels(fpop[,name])==2)
  }
  f2pop <- fpop[,which(f2vars)]
  
  
  manipulate(
    n1=slider(1,max.samp.sizes,initial=1,label="Sample Size n1"),
    n2=slider(1,max.samp.sizes,initial=1,label="Sample Size n2"),
    numvar=picker(as.list(names(npop))),
    facvar=picker(as.list(names(f2pop))),
{popvals <- pop[,numvar]
 popmin <- min(popvals)
 popmax <- max(popvals)
 breaker <- f2pop[,facvar]
 
 #Will need these later on:
 results <- tapply(popvals,breaker,mean)
 diffmeans <- results[1]-results[2]
 
 results2 <- tapply(popvals,breaker,var)
 maxsddiffmeans <- sqrt(results2[1]/max.samp.sizes+results2[2]/max.samp.sizes)
 sddiffmeans <- sqrt(results2[1]/n1+results2[2]/n2)
 
 #determine upper y limit for histogram of samples x1bar-x2bar, may decide not to use this
 ymax <- 1.2/(sqrt(2*3.1416)*maxsddiffmeans)
 
 p1 <- lattice::histogram(~popvals|breaker,type="density",
                 main="Population Distributions",
                 xlab=numvar,
                 col="red",
                 panel = function(x, ...) {
                   lattice::panel.histogram(x, ...)
                 }          )
 
 twopops <- split(popvals,breaker)
 pop1 <- twopops[[1]]
 pop2 <- twopops[[2]]
 reps <- sim.reps
 sim.means <- numeric(reps)
 for (i in 1:reps)  {
   samp1 <- sample(pop1,n1,replace=FALSE)
   samp2 <- sample(pop2,n2,replace=FALSE)
   sim.means[i] <- mean(samp1)-mean(samp2)
 }
 
 p2 <- lattice::histogram(~sim.means,type="density",
                 main=paste(reps,"Simulations of x1bar-x2bar"),
                 xlab="x1bar-x2bar",
                 panel = function(x, ...) {
                   lattice::panel.histogram(x, ...)
                 })
 
 print(p1,split=c(1,1,1,2), more=TRUE)
 print(p2,split=c(1,2,1,2))
 
 #Print out some numerical results
 quantities <- c("mu1-mu2","expected value of x1bar-x2bar","stdev of x1bar-x2bar")
 
 theoretical <-c(diffmeans,diffmeans,sddiffmeans)
 in.simulation <- c(NA,mean(sim.means),sd(sim.means))
 frm <- data.frame(theoretical,in.simulation)
 rownames(frm) <- quantities
 
 print(frm)
}
  )
}  #end of SampDist2Means

if(getRversion() >= "2.15.1")  utils::globalVariables(c("picker", "numvar","facvar","n1","n2"))
