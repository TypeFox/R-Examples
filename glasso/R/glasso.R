glasso=function(s, rho, zero=NULL, thr=1.0e-4,maxit=1e4,  approx=FALSE, penalize.diagonal=TRUE,start=c("cold","warm"), w.init=NULL,wi.init=NULL, trace=FALSE){
n=nrow(s)
# on return, cflag=1 means the procedure did not converge, =0 means  it did

BIG=10e9

if(!is.matrix(rho) & length(rho)!=1 & length(rho)!=nrow(s))
   {stop("Wrong number of elements in rho")}
if(length(rho)==1){
   if(rho==0){ 
 warning("With rho=0, there may be convergence problems if the input matrix 
 is not of full rank")
}}

if(is.vector(rho)){ rho=matrix(sqrt(rho))%*%sqrt(rho)}
if(length(rho)==1){rho=matrix(rho,ncol=n,nrow=n)}

if(!is.null(zero)){
   if(!is.matrix(zero)){ zero=matrix(zero,nrow=TRUE)}
   for(k in 1:nrow(zero)){
      i=zero[k,1]
      j=zero[k,2]
    rho[i,j]=BIG
    rho[j,i]=BIG
  }}

start.type=match.arg(start)
if(start.type=="cold"){
  is=0
 ww=xx=matrix(0,nrow=n,ncol=n)
}
if(start.type=="warm"){
    is=1
   if(is.null(w.init) | is.null(wi.init)){
       stop("Warm start specified: w.init and wi.init must be non-null")
     }
      ww=w.init
      xx=wi.init
  }

itrace=1*trace
ipen=1*(penalize.diagonal)
ia=1*approx
mode(rho)="double"
mode(s)="double"
mode(ww)="double"
mode(xx)="double"
mode(n)="integer"
mode(maxit)="integer"
mode(ia)="integer"
mode(itrace)="integer"
mode(ipen)="integer"
mode(is)="integer"
mode(thr)="double"


junk<-.Fortran("glasso",
n,
s,
rho,
ia,
is,
itrace,
ipen,
thr,
maxit=maxit,
ww=ww,
xx=xx,
niter=integer(1),
del=double(1),
ierr=integer(1),
PACKAGE="glasso"
)

ww=matrix(junk$ww,ncol=n)
xx=matrix(junk$xx,ncol=n)

if(junk$ierr!=0){stop("memory allocation error")}

critfun=
function(Sigmahati,s, rho, penalize.diagonal=TRUE){
d=det(Sigmahati)
temp=Sigmahati
if(!penalize.diagonal){diag(temp)=0}
val= -log(d)+sum(diag(s%*%Sigmahati))+sum(abs(rho*temp))
return(val)
}
crit=critfun(xx,s,rho,penalize.diagonal)
return(list(w=ww,wi=xx,loglik=-(n/2)*crit,errflag=junk$ierr,approx=approx, del=junk$del, niter=junk$niter))
}
