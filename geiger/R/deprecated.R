.deprecate=function(prev, curr, ...){
	warning(paste(sQuote(prev), "is being deprecated: use", paste(sQuote(curr), "instead", collapse=" or "), sep=" "), ...)
}

.defunctify=function(prev, curr, ...){
	warning(paste(sQuote(prev), "is no longer available: use", paste(sQuote(curr), "instead", collapse=" or "), sep=" "), ...)
}


.namespace=function(package){
    x=loadedNamespaces()
    package%in%x
}

runMedusa=function(phy, richness, estimateExtinction = TRUE, modelLimit = 20, cutAtStem=TRUE, startR=0.05, startE=0.5, ...){
	.deprecate("runMedusa", "medusa")
    message("Refer to documentation (?medusa) for information on summarizing output")
	medusa(phy, richness=richness, model=ifelse(estimateExtinction, "bd", "yule"), cut=ifelse(cutAtStem, "stem", "node"), partitions=modelLimit, init=c(r=startR, epsilon=startE), ...)
}

prune.extinct.taxa=function(phy, tol= .Machine$double.eps^0.5){
    .deprecate("prune.extinct.taxa", "drop.extinct")
	drop.extinct(phy, tol)
}

prune.random.taxa=function(phy, n){
    .deprecate("prune.random.taxa", "drop.random")
	drop.random(phy, n)
}

name.check=function(phy, data){
    if(is.null(names(data))) stop("'data' must be given as a vector or matrix with names")
	.deprecate("name.check", "geiger:::.treedata")
	.treedata(phy, data)
}

BDsim=function(nStart, b, d, times){
	.deprecate("BDsim", "sim.bd")
	sim.bd(b=b, d=d, n0=nStart, times=times)
}


birthdeath.tree=function(b, d, time.stop=0, taxa.stop=0, seed=0, print.seed=FALSE, return.all.extinct=TRUE){
	.deprecate("birthdeath.tree", "sim.bdtree")
	crit=.check.stoppingcrit(time.stop, taxa.stop)
	sim.bdtree(b=b, d=d, stop=crit, n=taxa.stop, t=time.stop, seed=seed, extinct=return.all.extinct) 
}

tip.disparity=function(phy, data, disp=c("avg.sq", "avg.manhattan", 
"num.states")){
	.deprecate("tip.disparity", "disparity")
	disparity(phy=phy, data=data, index=disp)
}

ic.sigma=function(phy, data){
	.deprecate("ic.sigma", "ratematrix")
	ratematrix(phy=phy, dat=data)
}

rate.estimate=function(time=0, n=0, phy=NULL, epsilon = 0, missing = 0, crown=TRUE, kendall.moran=FALSE){
	.deprecate("rate.estimate", c("bd.ms", "bd.km"))
	if(kendall.moran){
		bd.km(phy=phy, time=time, n=n, missing=missing, crown=crown)
	} else {
		bd.ms(phy=phy, time=time, n=n, missing=missing, crown=crown, epsilon=epsilon)
	}
}

node.leaves=function(phy, node){
	.deprecate("node.leaves", "tips")
	tips(phy, node)

}

getAncStates=function(x, phy){
    .defunctify("getAncStates", "phytools:::fastAnc", immediate.=TRUE)
    stop()
#    if(.namespace("phytools")) {
#        td=treedata(phy, x, sort=TRUE)
#       if(ncol(td$data)>1) res=apply(td$data, 2, function(y) fastAnc(td$phy, y)) else res=fastAnc(td$phy, td$data[,1])
#            attr(res, "phylo")=td$phy
#           return(res)
#    } else {
#        stop("'phytools' is unavailable")
#    }
}

deltaTree = function(phy, delta, rescale=TRUE){
    .deprecate("deltaTree", "rescale.phylo")
    f=rescale.phylo(phy, model="delta")
    f(delta=delta, rescale=rescale) 
}

lambdaTree = function(phy, lambda){
    .deprecate("lambdaTree", "rescale.phylo")
    rescale.phylo(phy, "lambda", lambda=lambda)
}

kappaTree = function(phy, kappa){
    .deprecate("kappaTree", "rescale.phylo")
    rescale.phylo(phy, "kappa", kappa=kappa)
}

ouTree = function(phy, alpha){
    .deprecate("ouTree", "rescale.phylo")
    rescale.phylo(phy, "OU", alpha=alpha)
}

tworateTree = function(phy, breakPoint, endRate){
    .deprecate("tworateTree", "rescale.phylo")
    mt=max(heights(phy))
    bk=breakPoint/mt
    f=rescale.phylo(phy, "nrate")
    f(time=bk, rate=endRate, rescale=FALSE)
}

linearchangeTree = function(phy, endRate = NULL, slope=NULL){
    .deprecate("linearchangeTree", "rescale.phylo")
    flag="'endRate' or 'slope' must be supplied"
    if (is.null(slope) && is.null(endRate)) stop(flag)
    if(!is.null(slope) && !is.null(endRate)) stop(flag)
    
    rootdepth <- max(heights(phy))
    toslope=function(endRate, rootdepth){
        (endRate-1)/rootdepth
    }
    if (is.null(slope)) {
        slope = toslope(endRate, rootdepth)
    }
    
    rescale.phylo(phy, "trend", slope=slope)
}

exponentialchangeTree = function(phy, endRate=NULL, a=NULL){
    .deprecate("exponentialchangeTree", "rescale.phylo")
    
    flag="'endRate' or 'a' must be supplied"
    if (is.null(a) && is.null(endRate)) stop(flag)
    if(!is.null(a) && !is.null(endRate)) stop(flag)

    rootdepth <- max(heights(phy))
    if (is.null(a)) a <- log(endRate)/rootdepth
    
    rescale.phylo(phy, "EB", a=a)
}

speciationalTree = function(phy){
    .deprecate("speciationalTree", "rescale.phylo")
    rescale.phylo(phy, "kappa", kappa=0)

}

rescaleTree = function(phy, totalDepth){
    .deprecate("rescaleTree", "rescale.phylo")
    rescale.phylo(phy, "depth", depth=totalDepth)

}


