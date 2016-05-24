###########################################
# for use with INLA
###########################################
# dat is all positive and we want density to have positive support as well
# remember to transform the density!
get.density.boundary.corrected=function(dat) {
    # first log transform dat to real line
    if (min(dat)<0) return (density(dat))
    else {
        dat1=log(dat)
        tmp=density(dat1)
        list (x=exp(tmp$x), y=tmp$y/exp(tmp$x))
    }
}

# return the mean, sd, CI of the transformed variable
getMeanSd=function (marginal, f="identity") {
    
    # interpolations suggested by Havard: do it on the original scale
    logtao=log(marginal[,1]); p.logtao=marginal[,2]*marginal[,1]
    fun = splinefun(logtao, log(p.logtao)) 
    h=0.001
    x = seq(min(logtao),max(logtao),by=h) 
    pmf = exp(fun(x))*h
#    sum (pmf) # Prob
#    x = seq(min(logtao)-sd(logtao)/2,max(logtao)+sd(logtao)/2,by=h) 
#    pmf = exp(fun(x))*h
#    sum (pmf) # Prob
#    x = seq(min(logtao)-sd(logtao),max(logtao)+sd(logtao),by=h) 
#    pmf = exp(fun(x))*h
#    sum (pmf) # Prob
#    x = seq(min(logtao)-sd(logtao)*2,max(logtao)+sd(logtao)*2,by=h) 
#    pmf = exp(fun(x))*h
#    sum (pmf) # Prob
    
    lower.boundary = rle(cumsum(pmf)>.025)$lengths[1]+1
    upper.boundary = rle(cumsum(pmf)>.975)$lengths[1]+1
#    if (pmf[lower.boundary]>.04) {
#        #stop ("getMeanSd(): pmf too large at lower boundary: "%+%pmf[lower.boundary])
#        return (rep(NA, 4))
#    }
#    if (pmf[upper.boundary]>.04) {
#        stop ("getMeanSd(): pmf too large at upper boundary"%+%pmf[upper.boundary])
#        return (rep(NA, 4))
#    }
    
    if (f=="identity") {
        func=function(x) { exp(x) }
    } else if (f=="inverse") {
        func=function(x) { exp(x)**-1 }
    } else if (f=="inversesqrt") {
        func=function(x) { exp(x)**-.5 }
    } else 
        stop ("getMeanSd(): function not supported "%+%f)
    
    mu = sum( pmf * func(x))
    stdev = (sum( pmf * func(x)**2 ) - mu**2) ** .5
    out = c("mean"=mu, "stddev"=stdev, 0,0 ) # we may run into problems with the boundary
    #out = c("mean"=mu, "stddev"=stdev, sort(func(x)[c(lower.boundary, upper.boundary)]) ); 
    names(out)[3]="2.5%"; names(out)[4]="97.5%";
    out
}

get.inla.f.param=function(arg) {
    if (arg$distr=="Gamma")
        inla.f.param=c(arg$shape, arg$rate)
    else if (arg$distr=="2DWishart")
        inla.f.param=c(arg$r, arg$R[1,1], arg$R[2,2], arg$R[1,2])
    else if (arg$distr=="3DWishart")
        inla.f.param=c(arg$r, arg$R[1,1], arg$R[2,2], arg$R[3,3], arg$R[1,2], arg$R[1,3], arg$R[2,3])
    else 
        stop ("get.inla.f.param(). distr not supported: "%+%arg$distr) 
    inla.f.param
}

getSamplesFileName=function (data.model,seed,chain,prior){
    seedstr=ifelse(is.null(seed) | is.na(seed),"","_seed"%+%seed)
    data.model%+%"/"%+%data.model%+%"_prior"%+%prior%+%"_chain"%+%chain%+%seedstr%+%".jags"
}

# plot marginal distribution from jags and inla
plot.marginals=function (.data.model, prior, seed=NA) {
    # first combine samples
    samples = NULL
    for (i in 1:.data.model@gibbs.chains) samples=rbind(.data.model@gibbs.samples[[i]])
    
    seedstr=""
    if (.data.model@is.sim) seedstr="_seed"%+%seed
    mypostscript(file=.data.model@name%+%"_prior"%+%prior%+%seedstr%+%".eps" )
    for (i in 1:length(.data.model@gibbs.watch.list)) {
        # this is a little precarious, but hey
        if (startsWith(.data.model@gibbs.watch.list[i],"sigma")) {
            jags.marginal = get.density.boundary.corrected ( 1/samples[,i]**2 ) 
        } else if (startsWith(.data.model@gibbs.watch.list[i],"rho")) {
            jags.marginal = get.density.boundary.corrected ( samples[,i] ) 
        } else if (startsWith(.data.model@gibbs.watch.list[i],"log")) {
            jags.marginal = get.density.boundary.corrected ( samples[,i] ) 
        } else if (.data.model@gibbs.watch.list[i]=="z") {
            jags.marginal = get.density.boundary.corrected ( samples[,i] ) 
        } else {
            next; # only do variance component parameters
        }
        
        hyper.marginal=.data.model@hyper.fit$marginals[[i]]
        # sometimes it is necessary to do transformation
        if (startsWith(.data.model@gibbs.watch.list[i],"log")) {
            hyper.marginal = cbind(log(hyper.marginal[,1]), hyper.marginal[,2]*hyper.marginal[,1])
        } else if (.data.model@gibbs.watch.list[i]=="z") {
            hyper.marginal = cbind(logit(.5+hyper.marginal[,1]/2), hyper.marginal[,2]*(.5-hyper.marginal[,1]**2/2))
        } else {
            hyper.marginal = hyper.marginal 
        }
        
        plot(jags.marginal$x, jags.marginal$y, type="l", xlab="", ylab="density", ylim=c(0,max(jags.marginal$y,hyper.marginal[,2])))
            #, xlim=range(c(marginal2[,1])) )
        lines(hyper.marginal[,1], hyper.marginal[,2], col=3)
#        plot(hyper.marginal[,1], hyper.marginal[,2], col=3, type="l", xlab="", ylab="density")
#        lines(jags.marginal$x, jags.marginal$y )
        legend(legend=c("mcmc","hyperpar"), col=c(1,3),x="topright", inset=0, bty="o", lty=1, cex=.6)
        title(main=.data.model@gibbs.watch.list[i])
    }
    dev.off()
}

load.jags.samples=function (.data.model, jags.samples.tmp, prior,seed=NA) {
    for (i in 1:.data.model@gibbs.chains) {
        load(getSamplesFileName(.data.model@name,seed,i,prior))
        .data.model@gibbs.samples[[i]] = jags.samples.tmp[,.data.model@gibbs.watch.list,drop=FALSE]
        print(getVarComponent(.data.model@gibbs.samples[[i]]))
    }
    .data.model
}
