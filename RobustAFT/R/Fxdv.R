"Fxdv" <-
function(v0,rr,n,Beta,tol,maxit){
# Fixed point algorithm for scale
args <- list(rr=rr,n=n,Beta=Beta)
v    <- v0; nit  <- 1
repeat{
  vo <- v
  v  <- AveS21w(vo,args)*vo; d <- v-vo
  if (nit==maxit | abs(d)<tol) break
  nit <- nit+1}
list(v=v,nit=nit,delta=d)}

