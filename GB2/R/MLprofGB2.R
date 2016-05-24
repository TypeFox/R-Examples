
# Maximum likelihood based on the profile log-likelihood

profml.gb2 <- function(z, w=rep(1, length(z)), method=1, hess = FALSE){

fnp <- function(x, z, w){
	a <- x[1]
	b <- x[2]
return(-proflogl.gb2(z, a, b, w))
}

grp <- function(x, z, w){
	a <- x[1]
	b <- x[2]
return(-profscores.gb2(z, a, b, w))
}

# Initial values of a and b under Fisk
x0 <- fisk(z, w)[1:2]  

opt1 <- optim(x0, fnp, grp, z, w, method="BFGS", control=list(parscale=x0,pgtol=1e-16), hessian=hess)
if (method != 2) return(list(opt1=opt1))
if (method == 2){
opt2 <- optim(x0, fnp, grp, z, w, method="L-BFGS-B", lower=0, control=list(parscale=x0,pgtol=0), hessian=hess)
return(list(opt2=opt2))
}
}