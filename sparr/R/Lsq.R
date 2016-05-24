Lsq <- function(x,uh,WIN){
    res <- KBivN_OPT(x=x,WIN=WIN,mu=c(uh[1],uh[2]),Sigma=matrix(c(uh[3],0,0,uh[3]),2,2),sq=F)
    return((2*res$val - res$u1^2*res$val - res$u2^2*res$val)^2)
}
