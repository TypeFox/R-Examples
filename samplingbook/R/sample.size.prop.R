sample.size.prop <-
function(e, P=0.5, N=Inf, level=0.95){    
# e = precision, wanted width of confidence interval
# p = expected proportion of events
# N = population size
# level = coverage probability for confidence intervals
            
  # input tests
  if(level<0 | level>1) stop("Wrong input: ", sQuote("level")," has to be probability between 0 and 1.")
  if(!is.numeric(N)) stop("Wrong input: ", sQuote("N")," is not a number or ", sQuote("Inf"))
  if(e<0) stop("Wrong input: precision ", sQuote("e")," has to be positive number.")
  if(e>=0.5) warning("Check input: precision bigger than 50% might be not useful for proportion confidence interval.")
  if(P<0 | P>1) stop("Wrong input: ", sQuote("P")," has to be probability between 0 and 1.")
  if(e==0 & N==Inf) stop("For approximate calculation (", sQuote("N")," much larger than ", sQuote("n"),") it must be ", sQuote("e > 0"),".")
  if((P<0.1 | P>0.9) & N<Inf){
    N <- Inf
    warning("Proportion ",sQuote("P")," is close to border. Using ", sQuote("N = Inf")," for calculations.\n")
  }
  if(e > P | e > (1-P)) stop("Wrong input: Estimate useless, if ", sQuote("e > P")," or ", sQuote("e > (1-P)"))
  q= qnorm((1+level)/2)
  # N large enough
  if(N==Inf) {
    n = q^2 * P*(1-P) / e^2
    n = ceiling(n)
  }  
  # For population size N
  else {
    n= P*(1-P) / (e^2/q^2 + P*(1-P)/N)
    n= ceiling(n)
  }
  # return argument
  ret <- list()
  ret$call <- list(e=e,P=P,N=N,level=level)
  ret$n <- n
  structure(ret,class="sample.size.prop")
}
