# Calcualtion of the vl, l=1,...,L-1 as a function of pl, l=1, ..., L
vofp.cgb2 <- function(pl){
ncomp <- length(pl)
vl <- rep(0, ncomp-1)
pL <- pl[ncomp]
vl <- log(pl[-ncomp]/pL)
return(vl)
}

# Calculation of the pl, l=1,...,L as a function of vl, l=1, ..., L-1
pofv.cgb2 <- function(vl){
pl <- exp(vl)/(1+sum(exp(vl)))
pL <- 1/(1+sum(exp(vl)))
pl <- c(pl,pL)
return(pl)
}

# Log-likelihood function 
logl.cgb2 <- function(fac, pl, w=rep(1, dim(fac)[1])){
    sw <- sum(w)
	mixt <- fac%*%pl
	logcomp <- log(mixt)
	logL <- sum(w*logcomp)/sw
	 return(logL)
}

# Score functions 
# Weighted mean of the scores
scores.cgb2 <- function(fac, pl, w=rep(1, dim(fac)[1])){
    sw <- sum(w)
    L <-length(pl)
	denom <- fac%*%pl
	num <- fac[,-L]
	midt <- num/as.vector(denom) - 1
	 if (L>2)  dlogL <- pl[-L]*colSums(midt*w)       
	 if (L==2) dlogL <- pl[-L]*sum(midt*w)         
	 return(dlogL/sw) 
}

# Maximum pseudo-likelihood estimation, GB2 as a compound distribution
ml.cgb2 <- function (fac, pl0, w = rep(1, dim(fac)[1]), maxiter = 100, fnscale=length(w)) 
{
    vl0 <- vofp.cgb2(pl0)
    fn <- function(vl, fac, w) {
        pl <- pofv.cgb2(vl)
        return(-logl.cgb2(fac, pl, w))
    }
    gr <- function(vl, fac, w) {
        pl <- pofv.cgb2(vl)
        return(-scores.cgb2(fac, pl, w))
    }
    opt <- optim(vl0, fn, gr, fac, w, method = "BFGS", control = list(maxit = maxiter, 
        fnscale = fnscale), hessian = FALSE)
    vlf <- opt$par
    plf <- pofv.cgb2(vlf)
    return(list(plf, opt))
}