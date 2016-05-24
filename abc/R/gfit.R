######################################################################
#
# gfit.R
#
# copyright (c) 2014-07-10, Katalin Csillery, Louisiane Lemaire, 
# Olivier Francois and Michael GB Blum
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Part of the R/abc package
# Contains: gfit, is.gfit, plot.gfit, summary.gfit
######################################################################


gfit=function(target, sumstat, nb.replicate, tol=.01, statistic=mean, subset=NULL, trace=FALSE){


    #checks : sumstat, nb.replicate and target
    if (missing(sumstat)){
        stop("'sumstat' must be supplied.", call.=FALSE)
    }   
    if (is.vector(sumstat)){sumstat=as.data.frame(sumstat)}
    if (missing(nb.replicate)){
        stop("'nb.replicate' must be supplied.", call.=FALSE)
    }   
    if (is.null(target)){
        stop("'target' must be supplied ", call.=FALSE)
    }

    #when subset is missing
    if (is.null(subset)){
    subset=rep(TRUE, dim(sumstat)[1])
    }


    #param
    param=data.frame(1:dim(sumstat)[1])


    #distance for the observed data
    res.abc.obs=abc(target=target, param=param, sumstat=sumstat, tol=tol, method="rejection", subset=subset)
    dist.obs=statistic(res.abc.obs$dist)


    #distances for the simulated data
    dist=function(obs){
        res=abc(target=sumstat[obs,], param=as.data.frame(param[-obs,]), sumstat=sumstat[-obs,], tol=tol, method="rejection", subset=subset[-obs])
        if (trace==TRUE){print(which(myseq==obs))}
        d=statistic(res$dist)
        return (d)
    }
    myseq=sample(c(1:dim(sumstat)[1]), replace=F, size=nb.replicate)
    dist.sim=sapply(myseq, dist)

    
    #return 
    gfit.out=list(dist.sim=dist.sim, dist.obs=dist.obs)
    class(gfit.out)="gfit"
    invisible(gfit.out)
}



is.gfit=function(x){
    if (inherits(x, "gfit")){
        return(TRUE)
    }else{ 
        return(FALSE)}
}



plot.gfit=function(x, breaks="Freedman-Diaconis", main="Histogramme of the null distribution", ...){

    if (!inherits(x, "gfit")) 
        stop("Use only with objects of class \"gfit.\"", call.=F)

    dist.sim=x$dist.sim
    dist.obs=x$dist.obs
    et=max(dist.sim)-min(dist.sim)
    inf=min(min(dist.sim), dist.obs)-0.1*et
    sup=max(max(dist.sim), dist.obs)+0.1*et
    
    hist(dist.sim, prob=T, xlab="Distance", xlim=c(inf, sup), breaks=breaks, main=main, ...)
    abline(v=dist.obs, col=4, pch=5, lwd=2)
    legend("topright", legend="Observed distance", cex=0.8, col=4, lty=1, lwd=2)
}



summary.gfit=function(object, ...){

    if (!inherits(object, "gfit")) 
        stop("Use only with objects of class \"gfit\".", call.=F)

    dist.sim=object$dist.sim
    dist.obs=object$dist.obs
    atypique=dist.sim[dist.sim>=dist.obs]
    pvalue=round(length(atypique)/length(dist.sim), 4)
    s.dist.sim=summary(dist.sim,...)

    return(list(pvalue=pvalue, s.dist.sim=s.dist.sim, dist.obs=dist.obs))
}

