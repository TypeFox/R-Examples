##
##  Wald Wolfowitz Runs Test
##
runs.test <- function(x, alternative="two.sided", threshold=median(x), pvalue="normal", plot=FALSE){
  # Performs the Runs Test for Randomness.
  #
  # Args:
  #   x: a numeric vector containing the data.
  #   alternative hypothesis, must be one of "two.sided" (default), "left.sided" or "right.sided"
  #   threshold: 
  #
  # Returns:
  #   statistic: the (normalized) value of the statistic test.
  #   n: the sample size, after the remotion of consecutive duplicate values.
  #   p.value: the asymptotic p-value.
  #
  dname <- deparse(substitute(x))
  if (alternative == "t"){alternative <- "two.sided"} 
  if (alternative == "l"){alternative <- "left.sided"}
  if (alternative == "r"){alternative <- "right.sided"}    
  if (alternative != "two.sided" & alternative != "left.sided" & alternative != "right.sided")
    {stop("must give a valid alternative")}
  # Remove NAs
  x <- na.omit(x)
  stopifnot(is.numeric(x))
  # Remove values equal to the level
  x <- x[x!=threshold]
  s <- sign(x-threshold)
  n1 <- length(s[s>0]) 
  n2 <- length(s[s<0])
  runs <- rle(s)
  r1 <- length(runs$lengths[runs$values==1])
  r2 <- length(runs$lengths[runs$values==-1])  
  n <- n1+n2
  mu <- 1 + 2*n1*n2/(n1+n2)
  vr <- 2*n1*n2*(2*n1*n2-n1-n2)/(n^2*(n-1))
  rr <- r1+r2
  #
  # Plot the data if requested by the user
  if (plot){
    plot((1:n)[s>0],x[s>0], xlim=c(1,n), ylim=c(min(x),max(x)), xlab="", ylab=dname)
    points((1:n)[s<0],x[s<0], col="red")
    abline(h=threshold, col=gray(.4))
    for (i in 1:(n-1)){
      if (s[i]*s[i+1]<0){abline(v=i+0.5, lty=2)}
    }
  }
  #
  # Computes the p-value
  pv <- 0
  if (pvalue == "exact"){    
    if (alternative=="two.sided"){
      pv1<-sum(druns(1:rr,n1,n2))
      pv2<-sum(druns(rr:(n1+n2),n1,n2))
      pv <- 2*min(pv1,pv2)
    }
    if (alternative=="left.sided"){pv<-sum(druns(2:rr,n1,n2))}
    if (alternative=="right.sided") {pv<-sum(druns(rr:(n1+n2),n1,n2))}    
  }
  if (pvalue=="normal"){
    pv0 <- pnorm((rr - mu) / sqrt(vr))
    if (alternative=="two.sided"){pv <- 2*min(pv0,1-pv0)}
    if (alternative=="left.sided"){pv <- pv0}
    if (alternative=="right.sided") {pv <- 1-pv0}
  }  
  if (alternative=="two.sided"){alternative<-"nonrandomness"}
  if (alternative=="left.sided"){alternative<-"trend"}
  if (alternative=="right.sided") {alternative<-"first-order negative autocorrelation"}
  #
  rval <- list(statistic = c(statistic=(rr - mu) / sqrt(vr)), p.value = pv, runs=rr, mu=mu, var=vr,  
               method = "Runs Test", data.name = dname, parameter=c(runs=rr, n1=n1,n2=n2,n=n), alternative=alternative)  
  class(rval) <- "htest"
  return(rval)
  
}  

# x<-c(68, 71, 69, 71, 70, 65, 63, 64, 65, 64, 67, 68, 66, 68, 66, 70)