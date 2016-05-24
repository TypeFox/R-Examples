#' @title Skewness in the Binomial Family of Distributions

#' @description An app to investigate how skewness in a binomial distribtution vanishes
#' when np is large enough.  Sample size is set at n = 50, but the user can vary p with a slider.
#' 
#' @rdname BinomSkew
#' @usage BinomSkew()
#' @return no value.  Graphical side-effects only.
#' @export
#' @author Homer White (hwhite0@@georgetowncollege.edu)
#' @examples
#' \dontrun{
#' if (require(manipulate)) BinomSkew()
#' }
BinomSkew <- function(){
  #Normal Approximation to Binomial (need n*p big enough)
  #Number of trials fixed.  You can vary p.  x-axis scale stays same.
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  manipulate(
    p=slider(0.01,0.99,step=0.01,initial=0.01,label="Success Chance p"),
    norm=checkbox(FALSE,"Show Normal Curve"),
{
  n <- 50 #fixed number of trials
  if (norm==TRUE){
    lower.lim <- min(0,qnorm(0.0001,mean=n*p,sd=sqrt(n*p*(1-p))))
    upper.lim <- max(n,qnorm(0.9999,mean=n*p,sd=sqrt(n*p*(1-p))))
  } else {
    lower.lim <- 0
    upper.lim <- n
  }
  plot(0:n,dbinom(0:n,size=n,prob=p),type="h",xlim=c(lower.lim,upper.lim),xlab="Number of Successes x",
       ylab="p(x)",main=paste("Binom(",n,", p):  Exp Successes = ",n*p,", Exp Failures = ",n*(1-p)))
  
  if (norm==TRUE) curve(dnorm(x,mean=n*p,sd=sqrt(n*p*(1-p))),n=100001,lwd=3,
                        from=lower.lim,to=upper.lim,col="red",add=TRUE)
  
}
  )
}
