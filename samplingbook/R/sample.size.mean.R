sample.size.mean <-
function(e, S, N=Inf, level= 0.95){

# S = standard deviation in population
# e = precision, wanted width of confidence interval
# N = population size
# level = coverage probability for confidence intervals

  # input tests
    if(level<0 | level>1) stop("Wrong input: ", sQuote("level")," has to be probability between 0 and 1.")
    if(e<0) stop("Wrong input: precision ", sQuote("e")," has to be positive number.")
    if(S<0) stop("Wrong input: standard deviation ", sQuote("S")," has to be positive number.")
    if(!is.numeric(N)) stop("Wrong input: ", sQuote("N")," is not a number or ", sQuote("Inf"),".")
    if(e==0 & N==Inf) stop("For approximate calculation (", sQuote("N")," much larger than ", sQuote("n"),") it must be ", sQuote("e > 0"))

  q= qnorm((1+level)/2)

  # N much larger than n
  if(N==Inf) { 
    n = q^2 * S^2 / e^2
    n = ceiling(n)
  }  
  # For population size N
  else  {
    n = S^2 / (e^2/q^2 + S^2/N)
    n = ceiling(n)
  } 
  # return argument
  ret <- list()
  ret$call <- list(e=e,S=S,N=N,level=level)
  ret$n <- n
  structure(ret,class="sample.size.mean")
}
