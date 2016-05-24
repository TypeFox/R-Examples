local(envir=.PBSmodEnv,expr={
locale = sys.frame(sys.nframe() - 1) # local environment

# ***********************************************************
# fib.C:
#   Compute Fibonacci numbers iteratively using a .C call to
#   C code
# Arguments:
#   n   - final nth fibonacci number to calculate
#   len - length of output vector with previous fibonacci numbers
# -----------------------------------------------------------
fib.C <- function(n=defaultN, len=defaultLen) {
  if (n<0) return(NA)
  if (len>n) len <- n
  retArr <- numeric(len)
  out <- .C("fibonacci", as.integer(n), as.integer(len),
            as.numeric(retArr))
  x <- out[[3]]
  return(x) }

# ***********************************************************
# fib.R:
#   A native R version of fib.C, used for comparison
# Arguments:
#   n   - final nth fibonacci number to calculate
#   len - length of output vector with previous fibonacci numbers
# ----------------------------------------------------------- 
fib.R=function(n=defaultN, len=defaultLen){
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
fib.init=function(){
  defaultN <- 200; tput(defaultN)
  defaultLen <- 10; tput(defaultLen)
}
#fib.init()
#print(fib.R(13))
#print(fib.C(13))

}) # end local scope
	