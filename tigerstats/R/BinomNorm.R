#' @title Binomial Distributions With Normal Approximation

#' @description An app to investigate the binomial family.
#' 
#' @rdname BinomNorm
#' @usage BinomNorm()
#' @return no value.  Graphical side-effects only.
#' @export
#' @author Homer White (hwhite0@@georgetowncollege.edu)
#' @examples
#' \dontrun{
#' if (require(manipulate)) BinomNorm()
#' }
BinomNorm <- function()  {
  #Normal approx to binomial.  Can vary n, but x-axis scale will change.
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  manipulate(
    p=slider(0.01,0.99,step=0.01,initial=0.20,label="Success Chance p"),
    n=slider(1,1000,initial=5,label="Number of Trials n"),
    norm=checkbox(FALSE,"Show Normal Curve"),
{
  lower.lim <- min(qbinom(0.0001,size=n,prob=p),qnorm(0.0001,mean=n*p,sd=sqrt(n*p*(1-p))))
  upper.lim <- max(qbinom(0.9999,size=n,prob=p),qnorm(0.9999,mean=n*p,sd=sqrt(n*p*(1-p))))
  max.binom <- max(dbinom(0:n,size=n,prob=p))
  max.norm <- max(dnorm(seq(lower.lim,upper.lim,by=1/101),mean=n*p,sd=sqrt(n*p*(1-p))))
  max.ylim <- max(max.binom,max.norm)
  ylims <- c(0,max.ylim)
  if (norm==TRUE){
    plot(0:n,dbinom(0:n,size=n,prob=p),type="h",xlim=c(lower.lim,upper.lim),xlab="Number of Successes x",
         ylim=ylims,ylab="p(x)",main=paste("Binom(",n,",",p,")"))
    curve(dnorm(x,mean=n*p,sd=sqrt(n*p*(1-p))),n=10001,lwd=3,
          from=lower.lim,to=upper.lim,col="red",add=TRUE)
  } else {
    plot(0:n,dbinom(0:n,size=n,prob=p),type="h",xlim=c(lower.lim,upper.lim),
         ylim=ylims,
         xlab="Number of Successes x",
         ylab="p(x)",main=paste("Binom(",n,",",p,")"))
  }
  
}
  )
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("n","p","x"))
