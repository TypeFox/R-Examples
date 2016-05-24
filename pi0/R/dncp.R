dncp=function(obj, ...) UseMethod("dncp")

dncp.parncpt=function(obj, fold=FALSE, ...)
{
    ans=function(x)dnorm(x,obj$mu.ncp,obj$sd.ncp)
    if(fold) fold(ans) else ans
}

dncp.parncpt2=function(obj, fold=FALSE, ...)
{
    ans=function(x) obj$tau*dnorm(x, obj$mu1.ncp, obj$sd1.ncp)+(1-obj$tau)*dnorm(x, obj$mu2.ncp, obj$sd2.ncp)
    if(fold) fold(ans) else ans
}

dncp.nparncpt=function(obj, fold=FALSE, ...) {
    d.ncp = function(xx) {
        xx = outer(xx, obj$all.mus, "-")
        xx = sweep(xx, 2, obj$all.sigs, "/")
        d = sweep(dnorm(xx), 2, obj$all.sigs, "/")
        drop(d %*% obj$beta)
    }
    if(fold) fold(d.ncp) else d.ncp    
}

dncp.sparncpt=function(obj, fold=FALSE, ...) 
{
    d.ncp=function(x) obj$par*dncp(obj$parfit)(x) + (1-obj$par)*dncp(obj$nparfit)(x)
    if(fold) fold(d.ncp) else d.ncp
}

dncp.nparncpp=function(obj, reflect=TRUE,...)
{
  .NotYetImplemented()
}


dncp.parncpF=function(obj,...)
{
    ans=function(x)dchisq(x/obj$gamma2, obj$data$df1, obj$delta0/obj$gamma2)/obj$gamma2
    ans
}

dncp.nparncpF=function(obj,...)        # p in the paper
{   ## depends on mus, sigs
	mu2s=sort(unique(obj$all.mus^2))
	gam2=obj$gam2
    d.ncp=function(xx)        # p in the paper
    {   ## depends on mu2s, gam2, obj
		schisq.x=outer(xx/gam2, mu2s/gam2, dchisq, df=obj$data$df1)/gam2
		ans=drop(schisq.x %*% obj$beta)
		if(any(xx==0) && obj$data$df1<2){
			ans[xx==0]=Inf
		}
		ans
    }
    d.ncp
}

dncp.sparncpF=function(obj,...) 
{
    d.ncp=function(x) obj$par*dncp(obj$parfit)(x) + (1-obj$par)*dncp(obj$nparfit)(x)
    d.ncp
}
