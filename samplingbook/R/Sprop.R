Sprop <-
function(y, m, n=length(y), N=Inf, level=0.95){

# y = binary data vector
# m = number of positive events
# n = sample size
# N = population size
# level = coverage probability for confidence intervals

### input treatment
    if(level<0 | level>1) stop("Wrong input: ", sQuote("level")," has to be probability between 0 and 1.")
    if(!is.numeric(N)) stop("Wrong input: ", sQuote("N")," is not a number or ", sQuote("Inf"))
    if(N>100000 & N!=Inf){
      N<-Inf
      warning("For ", sQuote("N > 100000")," using ", sQuote("N = Inf")," for calculation.")
    }
    if(n<=0) stop("Wrong input: ", sQuote("n")," has to be positive integer.")
    if(n > N) stop("Wrong input: sample size ", sQuote("n")," has to be smaller than population ", sQuote("N > 100000"))
    if(missing(m)){
     if(missing(y)) stop("Wrong input: ", sQuote("y")," or ", sQuote("m")," and ", sQuote("n")," has to be given.")
     else{
      if(any(is.na(y))) stop("Wrong input: No missing values in ", sQuote("y")," allowed.")
      if(!all(y==0|y==1)) stop("Wrong input: ", sQuote("y")," has to be a binary vector with values 0 and 1.")
      m <- sum(y)
     }
    }
    else{
     if(m<0) stop("Wrong input: ", sQuote("m")," has to be positive integer.")
     if((m - floor(m))>0) stop("Wrong input: ", sQuote("m")," has to be positive integer.")
     if(m > n) stop("Wrong input: positive events ", sQuote("m")," has to be smaller than sample size ", sQuote("n"),".")
     if(missing(y)) y <- NULL  
     else{
      if(m!=sum(y)) stop("Wrong input: ", sQuote("m")," is not the total number of events in ", sQuote("y"),)
     }
    }
### general calculations
  p <- m/n
  q <- qnorm((1+level)/2)
### finite populations
  if(N < Inf){ 
  var <- p*(1-p)/n * (n/(n-1))*(N-n)/N
  # calculation of hypergeometric confidence interval with finite population correction
    ugrh <- p-q*sqrt(var)
    ogrh <- p+q*sqrt(var)
    hynv.anteil <- c(ugrh,ogrh)
    hynv.anzahl = c(ceiling(ugrh*N),floor(ogrh*N))
  #exact procedure
    #lower border: start with smallest possible value (number 1 in sample)
    ugrex <- m
    #checkup if condition is complied with
    while (phyper(m-1,ugrex,N-ugrex,n)>(level+1)/2) { ugrex=ugrex+1 } 
    ugrex= max(m,ugrex-1)                                                                                 
    #upper border: total number - number 0 in sample
    ogrex=N-(n-m)
    while (phyper(m,ogrex,N-ogrex,n) < (1-level)/2) { ogrex=ogrex-1 }
    ogrex = min(N-(n-m), ogrex+1)
    exact.anteil = c(ugrex/N, ogrex/N)
    exact.anzahl = c(ugrex, ogrex)
  }
### infinite populations
  else{
  var =p*(1-p)/n
  # calculation of binomial confidence interval
    ugr  = p-q*sqrt(var)
    ogr  = p+q*sqrt(var)
    asym.anteil = c(ugr,ogr)
  # calculation of exact binomial confidence interval based on clopper pearson (r-function binom.test)
    btest = binom.test(m,n,conf.level=level)
    ucp  = btest$conf.int[[1]]
    ocp  = btest$conf.int[[2]] 
    cp.anteil = c(ucp,ocp)
  # calculation of binomial confidence interval by Wilson based on Agresti and Coull (1998)
    uac = ( p + q^2/(2*n) - q*sqrt(var + q^2/(4*n^2)) ) / (1+q^2/n) 
    oac = ( p + q^2/(2*n) + q*sqrt(var + q^2/(4*n^2)) ) / (1+q^2/n) 
    ac.anteil = c(uac,oac)
  }
# calculation of se
  if(var>0) se <- sqrt(var)
  else{
    se <- NA
    warning("Standard error is ", sQuote("NA"),", because calculations for variance of mean has been not positive. Confidence intervals may not be valid.")
  }  
### return argument
  ret <- list()
  ret$call <- list(y=y, m=m, n=n, N=N, level=level)
  ret$p <- p
  ret$se <- se
  if(N < Inf){
   ret$ci <- list(approx=hynv.anteil,exact=exact.anteil)
   ret$nr <- list(approx=hynv.anzahl,exact=exact.anzahl)
  }
  else{
   ret$ci <- list(bin=asym.anteil,cp=cp.anteil,ac=ac.anteil)
  }
  structure(ret, class= "Sprop")
}
