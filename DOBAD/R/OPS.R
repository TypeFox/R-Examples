###Charles Doss
### This is code relating to linear birth-death processes and the orthogonal polynomial
#   solution (OPS) as found by Karlin and McGregor.

## Don't know how to compute Pi0 with K&G formulas when no immigration.

###stupid way, could be coded in vector notation, but ...
process.prob.one <- function(t,lambda,mu,nu=0,X0=1, Xt,
                             #eps.t=.Machine$double.eps^.75,
                             eps.t=1e-10,
                             eps.params=1e-10,n=-111){
  myppo <- function(myt){process.prob.one.singlearg(t=myt,lambda=lambda,
                                                    mu=mu,nu=nu,X0=X0,
                                                    Xt=Xt,eps.t=eps.t,
                                                    eps.params=eps.params);}
  apply(as.matrix(t),1, myppo)
}


#####should code this
## process.prob.one.CTMC_PO_many <- function(ctmc, lambda,mu,nu,
##                              eps=.Machine$double.eps^.75 , n=-111){
  
## }

######should finish coding so t can be vector without using apply...
## process.prob.one <- function(t,lambda,mu,nu,X0, Xt,
##                              eps=.Machine$double.eps^.75 , n=-111){
##   (t<eps)* (1*(X0==Xt) + 0)
##   (t<eps)* 
##   if (t < eps) {
##     if (X0==Xt) return(1)
##     else return(0)
##   }
##   else if (nu >0){
##     beta <- nu/lambda;
##     if (lambda<mu){Pij.F(i=X0,j=Xt,t=t,L=lambda,m=mu,beta=beta)}
##     else Pij.E(i=X0,j=Xt,t=t,L=mu,m=lambda,beta=beta)
##   }
##   else {
##     if (X0==0) {if (Xt==0) return(1) else return(0)}
##     else if (Xt==0) return(process.prob.one.fft(t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt))
##     else if (lambda<mu) Pij.C(i=X0-1,j=Xt-1,t=t,L=lambda,m=mu,beta=2)
##     else Pij.D(i=X0-1,j=Xt-1,t=t,L=mu,m=lambda,beta=2)
##   }
## }


