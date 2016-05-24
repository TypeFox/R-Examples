local(envir=.PBSmodEnv,expr={
locale = sys.frame(sys.nframe() - 1) # local environment

# ***********************************************************
# fib2.C:
#   Compute Fibonacci numbers iteratively using a .Call() 
#   call to C code
# Arguments:
#   n   - final nth fibonacci number to calculate
#   len - length of output vector with previous fibonacci numbers
# -----------------------------------------------------------
fib2.C <- function(n=defaultN, len=defaultLen)
{
	retArr <- numeric(len)
	out <- .Call("fib2", as.integer(n), as.integer(len))
	return(out)
}

# ***********************************************************
# fib2.R:
#   A native R version of fib2.C, used for comparison
# Arguments:
#   n   - final nth fibonacci number to calculate
#   len - length of output vector with previous fibonacci numbers
# -----------------------------------------------------------
fib2.R=function(n=defaultN, len=defaultLen){
  if (n<0) return(NA)
  if (len>n) len <- n
  
  retArr <- numeric(len)
  
  xa=0; xb=1; xn=-1;
  
  for(i in 0:n){
    if(i <= 1)
      xn=i
    else{
      xn=xa+xb
      xa=xb
      xb=xn
    }
    j=i-n+len
    if(j>=0)
      retArr[j]=xn
  }
  return(retArr)
}

#initialization for testing
fib2.init=function(){
  defaultN <- 200; tput(defaultN)
  defaultLen <- 10; tput(defaultLen)
}

}) # end local scope
