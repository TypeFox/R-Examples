"Nwtv" <-
function(v0,rr,n,Beta,tol,maxit){
# Newton algorithm for scale
args <- list(rr=rr,n=n,Beta=Beta)
v    <- v0; nit  <- 1
repeat{
  f  <- AveS20w(v,args)
  f1 <- AveS2Pw(v,args)
  if (!is.finite(f) | !is.finite(f1)) {d <- NA; v <- NA; break}
  d  <- -f/f1
  if (is.nan(d) | !is.finite(d)) {d <- NA; v <- NA; break}
  v  <- v+d
# cat(nit,round(c(v,d),5),"\n")
  if (v<=0) {d <- NA; v <- NA; break}
  if (nit==maxit | abs(d)<tol) break
  nit <- nit+1}
list(v=v,nit=nit,delta=d)}

