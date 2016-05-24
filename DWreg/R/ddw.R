ddw<-function(x, q=exp(-1),beta=1)
{
if(q>1 | q <0)
	stop('q must be between 0 and 1', call. = FALSE)
if(beta <=0)
	stop('beta must be positive', call. = FALSE)
is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
if(beta==1)
 res<-dgeom(x,1-q)
else{
if(is.wholenumber(x))
  res<-q^(x^(beta))-q^((x+1)^(beta))
else
 {
res<-0
print("Warning message: non-integer value in ddw")
}}
return(res)
}