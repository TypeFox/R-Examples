crossProdLasso <- function(xtx,xty,rho,thr=1e-4,maxit=100,beta.init=NULL){
  return(list(beta=lasso_one(xtx,xty,rho,thr,maxit,beta.init)$beta))
}

lasso_one=function(w,ss, rho, thr=1.0e-4, maxit=100,trace=F, beta.init=NULL){
# does lasso fit of a single response variable  on predictors,
# via coordinate descent. Data is inner products w and ss

n=length(ss)

if(length(rho)==1){rho=rep(rho,n)}

itrace=1*trace

if(is.null(beta.init)){ beta.init=rep(0,n)}
mode(rho)="single"
mode(ss)="single"
mode(w)="single"
mode(n)="integer"
mode(maxit)="integer"
mode(itrace)="integer"
mode(thr)="single"
mode(beta.init)="single"

junk<-.Fortran("lasso7", rho, n, as.matrix(w), ss, thr, xx=beta.init, PACKAGE="scout")

return(list(beta=junk$xx))
}

