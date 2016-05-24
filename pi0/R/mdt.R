marginal.dt=function(obj,  ...) UseMethod('marginal.dt') ## marginal density of t

marginal.dt.parncpt=function(obj,  ...)
{

    single.df.ans=function(x, df) obj$pi0*dt(x, df )+(1-obj$pi0)*dtn.mix(x, df,obj$mu.ncp, obj$sd.ncp, ...)
    df.unique=sort(unique(obj$data$df))
    if(length(df.unique)==1) return(function(x)single.df.ans(x,obj$data$df[1]))

    function(x){    # discrete mixure of many distinct df's
        dftab=table(obj$data$df)
        prop=dftab/sum(dftab)
        sums=numeric(length(x))
        for(i in seq(along=dftab) ) sums=sums+prop[i]*single.df.ans(x, dftab[i])
        sums
    }
}

marginal.dt.nparncpt=function(obj,  ...)
{

    single.df.ans=function(x,df)
    {   sums=obj$pi0*dt(x, df)
        for(k in 1:length(obj$beta)) {
            tmp=dtn.mix(x,df, obj$all.mus[k], obj$all.sigs[k],...)
            if(any(tmp<0) || any(is.na(tmp))) {
                warning("Noncentral density unreliable. I switched to exact density function")
                tmp=dtn.mix(x,df, obj$all.mus[k], obj$all.sigs[k], approximation='none')
            }
            sums=sums+(1-obj$pi0)*obj$beta[k]*tmp
        }
        sums
    }
    df.unique=sort(unique(obj$data$df))
    if(length(df.unique)==1) return(function(x)single.df.ans(x, obj$data$df[1]))

    function(x){    # discrete mixure of many distinct df's
        dftab=table(obj$data$df)
        prop=dftab/sum(dftab)
        sums=numeric(length(x))
        for(i in seq(along=dftab) ) sums=sums+prop[i]*single.df.ans(x, dftab[i])
        sums
    }
}

marginal.dt.sparncpt=function(obj, ...){ ## to be implemented
    ans=function(x) obj$par*marginal.dt(obj$parfit,...)(x) + (1-obj$par)*marginal.dt(obj$nparfit,...)(x)
    ans
}
