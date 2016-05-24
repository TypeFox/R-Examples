signmedian.test <-
function (x,mu=0,alternative = c("two.sided", "less", 
          "greater"), conf.level = 0.95,exact=TRUE, ...)
{
  DNAME <- deparse(substitute(x))
  if(!is.numeric(x) ||length(x)<2)
    stop("'x' must be numeric and a vector")
  if((!is.numeric(mu) || length(mu)>1 || is.na(mu) )&& !missing(mu))
    stop("'mu' must be numeric and a single number")
  
  alternative <- match.arg(alternative)
  if (!missing(mu) && ((length(mu) > 1L) || !is.finite(mu))) 
    stop("'mu' must be a single number")
  if (!((length(conf.level) == 1L) && is.finite(conf.level) && 
          (conf.level > 0) && (conf.level < 1))) 
    stop("'conf.level' must be a single number between 0 and 1")
  
  if (missing(exact)) 
    exact <- TRUE
  if (!exact) 
    exact <- FALSE
  
  mx <- median(x)
  N<-length(x)
  
  x0 <- x[x!=mu]
  B <- length(x0[x0>mu])
  n <- length(x0)
  
  if(exact) {
    b.t <- binom.test(B,n,alternative=alternative,conf.level=conf.level)
    PVAL <- b.t$p.value
    a <- (1-conf.level)/2
    p <- NULL
    for(i in 0:round(N/2)){
      p[i+1]=pbinom(i,n,0.5)
    }
    A <- rep(a,round(N/2)+1)
    C <- which.min(abs(p-A))
    p.L<-sort(x)[C]
    p.U <-sort(x)[N+1-C]
    conf.level<-1-2*p[C]
    CINT <- c(p.L,p.U)
  }
  else{
    E<-n/2;
    V<-n/4;
    bl<-(n-B-E+0.5)/sqrt(V);
    bg<-(B-E-0.5)/sqrt(V);
    
    PVAL <- switch(alternative, less = pnorm(bl), greater = 1-pnorm(bg),two.sided=ifelse((pnorm(bl)>=0.5),1,(2*pnorm(bl))))
    alpha<-1-conf.level
    C <- floor(N/4-pnorm(alpha/2)*sqrt(N/4))
    p.L <- sort(x)[C]
    p.U <- sort(x)[N+1-C]
    conf.level <- pnorm((N-1-C-N/2)/sqrt(N/4))-pnorm((C-N/2)/sqrt(N/4))
    CINT <- c(p.L,p.U)
    
  }
  nm_alternative <- switch(alternative, two.sided = "the median of x is not equal to mu", 
                           less = "the median of x is less than mu", greater = "the median of x is greater than mu")
  STATISTIC <- switch(alternative, two.sided = length(x[x!=mu]), 
                      greater = B, less = n-B)
  names(STATISTIC) <- switch(alternative, two.sided = paste0("#(x!=",mu,")"), 
                             greater = paste0("#(x>",mu,")"), less = paste0("#(x<",mu,")"))
  attr(CINT, "conf.level") <- conf.level
  ESTIMATE <- mx
  names(x) <- "x"
  names(mu) <- "mu"
  names(ESTIMATE) <- "point estimator"
  if(exact) {
    RVAL <- list(statistic = STATISTIC,parameter=mu, p.value = PVAL, alternative = nm_alternative, 
                 conf.int=CINT,estimate = ESTIMATE, method = "Exact sign test", data.name = DNAME)}
  else {RVAL <- list(statistic = STATISTIC,parameter=mu, p.value = PVAL, alternative = nm_alternative, 
                     conf.int=CINT,estimate = ESTIMATE, method = "Asymptotic sign test(with continuity correction)", data.name = DNAME)}
  class(RVAL) <- "htest"
  return(RVAL)
}
