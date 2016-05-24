loqo <- function(c, H, A, b, l, u, r, sigf=7, maxiter=10000, margin=0.05, bound=10, verb=0,restart = 0){

    n <- dim(H)[1]
    m <- dim(A)[1]

	primal <- double(3*n)
	dual   <- double(2*n+1)
	
	.C("R_pr_loqo",n=n,m=m,linear=c,hessian=H,constr.mat=A,constr.vec=b,lower=l,upper=u,
	primal=primal,dual=dual,verbosity=as.integer(verb),sigf=as.integer(sigf),maxit=maxiter,
	margin=margin,bound=bound,restart=as.integer(restart),convergence = as.integer(0))
}