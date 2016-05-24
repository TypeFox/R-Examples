inv.surv <-
function(base.dist, parms, alpha, val){
	uniroot(surv.dist, lower=0,upper=100000, base.dist=base.dist, parms=parms, alpha=alpha, xbeta=val[1], res=val[2])$root
}
