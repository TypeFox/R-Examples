#' @title Circular Systematic Sampling
#'
#' @description Circular systematic sampling.
#'
#' @param N The population size.
#' @param n The sample size.
#'
#' @keywords Sampling
#' @export
#' @examples
#' circularSampling(500, 30)
`circularSampling` <- function(N=500, n=30){
  k=round(N/n, 0)
  r = sample(1:N, 1)
  f=n-1
  circ.syst=r
  for (j in 1:f) {
    if (r + j*k <= N)
      circ.syst = c(circ.syst, r+j*k)
    else
      (circ.syst = c(circ.syst, r+j*k-N))
  }
  print(paste("Your Circular Systematic Sample is: ", return(circ.syst)))
}
NULL




#' @title Linear Systematic Sampling
#'
#' @description Linear systematic sampling.
#'
#' @param N The population size.
#' @param n The sample size.
#'
#' @keywords Sampling
#' @export
#' @examples
#' linearSampling(500, 30)
#'
`linearSampling` <-function(N=500, n=30){
  k=round(N/n, 0)
  r = sample(1:k, 1)
  syst.samp= seq(r, r+k*(n-1), k)
  print(paste("Your Linear Systematic Sample is:  ", return(syst.samp)))
}
NULL



# Suppose that a book of N = 500 pages will be examined in order to estimate the total number of errors. For our purpose, we will randomly select a sample of n = 30 pages.
# Explain how you would select the sample if you use srswor;
#Explain how you would select the sample if you use linear systematic
# sampling ;
# Explain how you would select the sample if you use circular systematic
# sampling.




#' @title Replicated Simple Random Sample (SRS)
#'
#' @description Replicated Simple Random Sample (SRS).
#'
#' @param N The population size.
#' @param n The sample size.
#' @param g Number of independent sub-samples, each containing m = n/g units. Notice that m has to be a multiple of n and g.
#'
#' @keywords Sampling
#' @export
#' @examples
#' replicatedSampling(500, 30, 6)
#'
`replicatedSampling` = function(N=500, n=30, g=6){
  m=n/g
  samp=sample(1:N, n);
  samp=data.frame(matrix(samp, nrow=m));
  colnames(samp) <- paste0("g", 1:g, sep="")
  return(samp)
}
NULL




#' @title Replicated-Systematic Random Sampling
#'
#' @description Replicated-Systematic Random Sampling.
#'
#' @param N The population size.
#' @param n The sample size.
#' @param g Number of independent sub-samples, each containing m = n/g units. Notice that m has to be a multiple of n and g.
#'
#' @examples
#' linearReplicatedSampling(500, 30, 6)
#' @export
`linearReplicatedSampling` <- function(N=500, n=30, g=6){
m=round(n/g,0)
k=round(N/m, 0)
syst.samp=NULL
for (i in 1:g){
r = sample(1:k, 1)
syst.samp = c(syst.samp, seq(r, k*m, k))
}
syst.samp=matrix(syst.samp, nrow=m)
print(paste("Your Replicated Linear-Systematic Sample is:  ", return(syst.samp)))
}
NULL





#' @title Replicated Circular-Systematic Sampling
#' @description Replicated circular systematic sampling.
#' @param N The population size.
#' @param n The sample size.
#' @param g Number of independent sub-samples, each containing m = n/g units. Notice that m has to be a multiple of n and g.
#'
#' @examples
#' circularReplicatedSampling(500, 30, 6)
#' @export
`circularReplicatedSampling`<- function(N=500, n=30, g=6){
m=round(n/g,0)
k=round(N/m, 0)
circ.syst=NULL
for (i in 1:g){
  r = sample(1:N, 1)
  f=m-1
  for (j in 0:f){
    if (r + j*k <= N)
      circ.syst = c(circ.syst, r+j*k)
    else
        (circ.syst = c(circ.syst, r+j*k-N))
    }
}
circ.syst=matrix(circ.syst, nrow=m)
print("Your Replicated Circular-Systematic Sample is:  ", return(circ.syst))
}
NULL







#' @title Simple Sample Size for Surveys
#' @description Compute sample size for surveys.
#' @param p The proportion.
#' @param delta The error size.
#' @param popsize An integer for the population size.
#' @param deff An intger for the deff.
#' @param alpha The level of alpha/significance.
#'
#' @export
#' @examples
#' # Comercial public opinion samples in Brazil:
#' sampleSize(p=.50, delta=.03)
#' sampleSize(p=.50, delta=.02)
`sampleSize` <- function(p, delta = "auto", popsize=NULL, deff=1, alpha = .05){
  q <- 1-p
  pq <- cbind(p, q)
  minpq <- apply(pq, 1, min)
  if(any(delta=="auto")){
    delta <- ifelse(minpq >= .3, 0.1, ifelse(minpq >= .1, .05, minpq/2))
  }
  if(any(p >= 1) | any(delta >= 1) | any(popsize < 2) )
    stop("Proportion and delta both must < 1. Popsize must be >=2")
  else {
    n1 <- stats::qnorm(1-alpha/2)^2*p*(1-p)/delta^2
    if (!is.null(popsize)){
      n1 = n1/(1+n1/popsize)
    }
    if (deff != 1) {
      n1 = n1*deff }
  }
  deff1 <- deff
  if(deff==1) deff1 <- NULL
  table1 <- cbind(p, popsize, deff1, delta, round(n1))
  colnames(table1)[colnames(table1)=="deff1"] <- "deff"
  colnames(table1)[ncol(table1)] <- "n"
  returns <- list(p = p, delta=delta, popsize=popsize, deff=deff,
                  alpha = alpha, n1=n1, minpq=minpq,
                  table = as.data.frame(table1))
  class(returns) <- c("sampleSize", "list")
  returns
}
NULL

### print sample.size
#' @export
`print.sampleSize` <- function(x, ...)
{
  if(nrow(x$table) < 6){
    cat("\n")
    cat("Sample size for survey.","\n")
    cat("Assumptions:", "\n")
    cat("  Proportion       =", x$p, "\n")
    cat("  Confidence limit =", round((1-x$alpha)*100), "%","\n")
    cat("  Delta            =", x$delta, "from the estimate.", "\n")
    if (!is.null(x$popsize)){
      cat("  Population size  =", x$popsize, "\n")
    }
    if (x$deff != 1) {
      cat("  Design effect    =", x$deff, "\n")
    }
    cat("\n")
    cat("  Sample size      =", round(x$n1), "\n")
    cat("\n")
  }else{
    cat("Sample size for surveys.","\n")
    cat("Assumptions:", "\n")
    if(length(x$alpha) == 1) cat("  Confidence limit =", round((1-x$alpha)*100), "%","\n")
    if(length(x$delta) == 1)	cat("  Delta            =", x$delta, "from the estimate.", "\n")
    print(x$table, rownames=FALSE)
  }
}
NULL







#' @title Calculate and plot power of a sample
#'
#' @description Calculates and plots power of a sample z-test of a sample mean mu1
#'  against a population mean \code{mu0} (H_{0}: mu0 = mu1, H_{1}: mu0 <> mu1).
#'
#' @param mu0 This should be the "known" mean value for your population.
#' @param mu1 This should be the "expected" mean value from your sample. The delta between mu(0) and mu(1) is what you should consider a significant difference for the test.
#' @param n The sample size.
#' @param sigma This should be the known sigma (standard deviation) for the population.
#' @param alpha  This is the significance level, default is alpha(twosided) = .05.
#'
#' @details
#' \code{sample.power} calculates the power of a one-sample z-test (twosided)
#' and plots the density distributions under the assumption of of H_{0}: m = mu0 and
#' H_{1}: m = mu1. The rejection regions of H_{0} (alpha) are colored blue, while the rejection region of H_{1} (beta) is colored red.
#'
#' @return
#' \code{n} the sample size;
#' \code{sigma} the standard deviation;
#' \code{SE} the standard error of the mean;
#' \code{mu0} the mean of H_{0} in the population;
#' \code{mu1} the sample mean;
#' \code{mean.crit} the critical value of sample mean to achieve significance;
#' \code{ES} the population "effect" size gamma;
#' \code{delta} the effect size delta (Cohen);
#' \code{alpha} the significance level alpha (twosided);
#' \code{power} the power (1-beta).
#'
#' @examples
#' samplePower(mu0=68, mu1=69, sigma=3.1, n=100)
#' ## gives a power of .90
#'
#' @export
`samplePower` <- function(mu0=0, mu1=0, sigma=1, n=100, alpha=.05) {
  gamma = (mu1-mu0)/sigma
  delta = gamma*sqrt(n)
  se.mean = sigma/sqrt(n)

  z0L = mu0-stats::qnorm(1-alpha/2)*se.mean
  z0U = mu0+stats::qnorm(1-alpha/2)*se.mean
  z0min = mu0-3.5*se.mean
  z0max = mu0+3.5*se.mean
  z1min = mu1-3.5*se.mean
  z1max = mu1+3.5*se.mean
  if (mu1 > mu0) mean.crit=z0U
  else mean.crit=z0L

  if (mu1 != mu0) {
    cat('\nOne-sample z-test power calculation\n\n')
    cat('          n =',n,'\n')
    cat('          \u03C3 =',round(sigma,3),'\n')
    cat('   SE(mean) =',round(se.mean,3),'\n')
    cat('         \u03BC0 =',round(mu0,3),'\n')
    cat('         \u03BC1 =',round(mu1,3),'\n')
    cat('  mean.crit =',round(mean.crit,3),'\n\n')
    cat('effect size =',round(gamma,3),'\n')
    cat('      delta =',round(delta,3),'\n')
    cat('  sig.level =',round(alpha,3),'\n')
    cat('      power =',round(1-stats::pnorm(stats::qnorm(1-alpha/2,mean=0,sd=se.mean),
                                             mean=abs(mu1-mu0),sd=se.mean),3),'\n')
    cat('alternative = two-sided\n\n')
    if (n < 30) cat('Warning: N too small for z test!\n\n')
  }

  # Normal curve H0:

  graphics::curve(stats::dnorm(x,mean=mu0,sd=se.mean),from=z0min,to=z0max,
                  xlim=range(z0min,z0max,z1min,z1max),ylab='density',col="blue",
                  xlab=paste('\u03B1 = ',round(alpha,3),' (two-sided)',
                             ', \u03B2 = ',round(stats::pnorm(stats::qnorm(1-alpha/2,mean=0, sd=se.mean),mean=abs(mu1-mu0), sd=se.mean),3),
                             ', power = ',round(1-stats::pnorm(stats::qnorm(1-alpha/2,mean=0, sd=se.mean),mean=abs(mu1-mu0), sd=se.mean),3),sep=""),
                  main=paste('Power of z-Test of the Mean of a Single Population\n',
                             'n = ',n,', \u03BC0 = ',mu0,', \u03BC1 = ',mu1,', \u03C3 = ',sigma, sep=""))
  # Acceptance region H0:
  x=seq(z0L,z0U,min(.001,1/n))
  graphics::polygon(c(z0L,x,z0U),c(0,stats::dnorm(x,mean=mu0,sd=se.mean),0),col="lightyellow")

  # Normal curve H1:
  graphics::curve(stats::dnorm(x,mean=mu1,sd=se.mean),from=z1min,to=z1max,add=TRUE,col="red")

  graphics::text(mu0,0, '\u03BC0',pos=1,offset=.15,cex=.8)
  graphics::text(mu1,0, '\u03BC1',pos=1,offset=.15,cex=.8)

  if (mu1 > mu0){
    graphics::text(z0U,0,round(z0U,2),pos=1,offset=.15,cex=.8)

    # Acceptance region H1:
    x=seq(z0U,z1max,min(.001,1/n))
    graphics::polygon(c(z0U,x,z1max),c(0,stats::dnorm(x,mean=mu1,sd=se.mean),0),
                      col="lightyellow")

    # Rejection region H0:
    x=seq(z0U,z0max,min(.001,1/n))
    graphics::polygon(c(z0U,x,z0max),c(0,stats::dnorm(x,mean=mu0,sd=se.mean),0),density=20,
                      col="blue")
    x=seq(z0min,z0L,min(.001,1/n))
    graphics::polygon(c(z0min,x,z0L),c(0,stats::dnorm(x,mean=mu0,sd=se.mean),0),density=20,
                      col="blue")

    # Rejection region H1:
    if (z1min < z0U) {
      x=seq(z1min,z0U,min(.001,1/n))
      graphics::polygon(c(z1min,x,z0U),c(0,stats::dnorm(x,mean=mu1,sd=se.mean),0),density=20,
                        angle=135,col="red")
    }
  }
  else if (mu1 < mu0) {
    graphics::text(z0L,0,round(z0L,2),pos=1,offset=.15,cex=.8)

    # Acceptance region H1:
    x=seq(z1min,z0L,min(.001,1/n))
    graphics::polygon(c(z1min,x,z0L),c(0,stats::dnorm(x,mean=mu1,sd=se.mean),0),
                      col="lightyellow")

    # Rejection region H0:
    x=seq(z0U,z0max,min(.001,1/n))
    graphics::polygon(c(z0U,x,z0max),c(0,stats::dnorm(x,mean=mu0,sd=se.mean),0),density=20,
                      col="blue")
    x=seq(z0min,z0L,min(.001,1/n))
    graphics::polygon(c(z0min,x,z0L),c(0,stats::dnorm(x,mean=mu0,sd=se.mean),0),density=20,
                      col="blue")

    # Rejection region H1:
    if (z0L < z1max) {
      x=seq(z0L,z1max,min(.001,1/n))
      graphics::polygon(c(z0L,x,z1max),c(0,stats::dnorm(x,mean=mu1,sd=se.mean),0),density=20,
                        angle=135,col="red")
    }
  }
  else cat('Error: \u03BC1 must differ from \u03BC0!\n')
}
NULL
