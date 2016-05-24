#' @title Difference of Two Sample Proportions

#' @description An app to explore the sampling distribution of the difference of two sample proportions.
#' 
#' @rdname SampDist2Props
#' @usage SampDist2Props(form,data,max.sample.sizes=100,sim.reps=1000)
#' @param form An object of class formula, of the form ~x+y where x and y are factors supplied by:
#' @param data A dataframe, representing the imaginary population.  In the formula, both factors should have exactly two levels.
#' The variable x represents the explanatory variable.
#' @param max.sample.sizes Maximum sample sizes allowed on the sliders.
#' @param sim.reps  Number of samples to construct the empirical distribution. 
#' @return Graphical and numerical output.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' \dontrun{
#' data(imagpop)
#' SampDist2Props(~sex+cappun,data=imagpop)
#' }
SampDist2Props <-
function(form,data,max.sample.sizes=100,sim.reps=1000) {
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  
  #pop should be a data frame with numerical and factor variables
  
  #We are NOT gonna monkey around with missing values:
  data <- data[complete.cases(data),]
  
  prsd <- with(data,ParseFormula(form))
  
  Explanatory.Variable <- as.character(prsd$rhs)[2]
  Response.Variable <- as.character(prsd$rhs)[3]
  
  breaker <- data[,Explanatory.Variable]
  respvals <- data[,Response.Variable]
  
  twopops <- split(respvals,breaker)
  pop1 <- twopops[[1]]
  pop2 <- twopops[[2]]
  
  explanatory <- breaker
  response <- respvals
  results <- table(explanatory,response)
  print(results)
  cat("\n")
  
  N1 <- sum(results[1,])
  N2 <- sum(results[2,])
  p1 <- results[1,1]/N1
  p2 <- results[2,1]/N2
  
  if (max.sample.sizes>min(N1,N2))  {
    stop("Desired sample sizes must be less than pop sizes:  choose lower value for max.sample.sizes argument.")
  }
  
  
  
  manipulate(
    n1=slider(1,max.sample.sizes,initial=1,label="Sample Size n1"),
    n2=slider(1,max.sample.sizes,initial=1,label="Sample Size n2"),
    curvetype=picker("None","Density Estimate","Theoretical Normal"),
{#Will need these later on:
 
 big <- max.sample.sizes
 sddiffprops <- sqrt(p1*(1-p1)/n1*(N1-n1)/(N1-1)+p2*(1-p2)/n2*(N2-n2)/(N2-1))
 minsddiffprops <- sqrt(p1*(1-p1)/big*(N1-big)/(N1-1)+p2*(1-p2)/big*(N2-big)/(N2-1))
 #determine upper y limit for histogram of samples
 ymax <- 1.2/(sqrt(2*3.1416)*minsddiffprops)
  
 reps <- sim.reps
 sim.pdiffs <- numeric(reps)
 for (i in 1:reps)  {
   samp1 <- sample(pop1,n1,replace=FALSE)
   res.1 <- table(samp1)
   samp2 <- sample(pop2,n2,replace=FALSE)
   res.2 <- table(samp2)
   sim.pdiffs[i] <- res.1[1]/n1-res.2[1]/n2
 }
 
 hist(sim.pdiffs,freq=FALSE,xlim=c(-1,1),
      xlab=expression(hat(p)[1]-hat(p)[2]),
      ylim=c(0,ymax),
      col="blue",
      main=paste(reps,"simulations"))
 
 if(curvetype=="Density Estimate") {
   d <- density(sim.pdiffs)
   lines(d$x,d$y,col="red")
 }
 
 if (curvetype=="Theoretical Normal")  {
   curve(dnorm(x,mean=(p1-p2),sd=sddiffprops),
               xlim=c(-1,1),add=TRUE,col="red",n=1001)
 }
 
 #Print out some numerical results
 quantities <- c("p1-p2","expected value of p1hat-p2hat","stdev of p1hat-p2hat")
 
 theoretical <-c(p1-p2,p1-p2,sddiffprops)
 in.simulation <- c(NA,mean(sim.pdiffs),sd(sim.pdiffs))
 frm <- data.frame(theoretical,in.simulation)
 rownames(frm) <- quantities
 
 print(frm)
 cat("\n")
}
  )
  
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("picker", "curvetype","n1","n2"))
