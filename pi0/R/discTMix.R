discTMix= function(tstat, n1=10, n2=n1, nq, p0, p1, D, delta, paired=FALSE,
                    tbreak, ext=TRUE, threshold.delta=0.75, ...)
#
# Name: tMixture
# Desc: fit a mixture of t-components to a vector of t-statistics
# Auth: AP based on earlier tMixture4 170605
#
# Chng: 200605 AP added pairwise option
#       020905 AP changed pairwise to paired
#       021105 AP set default paired=FALSE
#       260106 AP changed np to nq, added threshold for components
#       060206 AP fixed eqDelta for nq=2, added check for nq<2
#                 remove dimension-attribute to be compatible with tstatistics()
#	 	07/04/15 LQ truncated output from ff; edited final output list
{
    # We want to deal with a vector of t-statistics, whch is unfortunately 
    # not what tstatistics() delivers
    if (!is.null(dim(tstat))) {
        if (ncol(tstat)==1) {
            tstat = tstat[,1] # silently
        } else if ("tstat" %in% colnames(tstat)) {
            tstat = tstat[,"tstat"] # silently
        } else {
            warning("only first column of tstat is used")
            tstat = tstat[,1]
        }
    }
    # Prepare essential and available parameters
    n = length(tstat)
    if (paired) {
        if (n1 != n2) {
            stop("paired requires equal sample sizes n1 and n2")
        }
        df = n1-1
        sf = sqrt(n1)
    } else {
        N = n1+n2
        df = N-2
        sf = sqrt(n1*n2/N)
    }        
    
    # Helper functions to make filling in other starting values more readable
    eqDelta = function(nq, tstat) {
                    nqhalf = floor((nq-1)/2)
                    delta  = if (nqhalf==0) NULL else c(-nqhalf:-1, 1:nqhalf)
                    if (nqhalf < (nq-1)/2) {
                        tt = table(tstat<0)
                        delta = c(delta, if (tt[1]<tt[2]) -0.5 else 0.5)
                    }
                    sort(delta)
            }
    eqP1 = function(nq, p0) rep((1-p0)/(nq-1), nq-1)
    # Find some kind of reasonable nq or break    
    if (missing(p1) & missing(D) & missing(delta)) {
        if (missing(nq)) {
            stop("Please use any reasonable combination of nq, p1, D, and delta",
                 " to specify an initial mixture model")
        }  else if (nq < 2) {
            stop("Specify at least two components")
        }
    } else {
        if (!missing(nq)) {
            warning("argument nq redundant - ignored")
        }
        if (!missing(p1)) {
            nq = length(p1) + 1
        } else if (!missing(delta)) {
            nq = length(delta)+1
        } else if (!missing(D)) {
            nq = length(D) + 1
        }
    }
    # Set the probabilities
    if (missing(p1)) {
        if (missing(p0)) {
            p0 = 1/nq
        }
        p1 = eqP1(nq, p0)
    } else {
        p1 = abs(p1)
        if (missing(p0)) {
            if (sum(p1)<1) {
                p0 = 1-sum(p1)
            } else {
                stop("Please specify a valid p0 either explicitly or as complement",
                     " of p1")
            }
        }
        # Renormalize
        pp = p0 + sum(p1); p0 = p0/pp; p1 = p1/pp
    }
    # Set the effects
    if (missing(D)) {
        if (missing(delta)) {
            delta = eqDelta(nq, tstat)
        }
        D = delta/sf
    } else {
        if (!missing(delta)) {
            warning("D is redundant when delta is given - D will be ignored")
            D = delta/sf
        } else {
            delta = D*sf
        }
    }
    # Check everybody agrees
    if (length(p1) != length(delta)) {
        stop("Unequal number of components for probabilities and effects")
    }
       
    # Tabulate
    if (missing(tbreak)) {
        tbreak = floor(sqrt(n))
    }
    if (length(tbreak)==1) {
        rr = range(tstat, na.rm=TRUE)
        tbreak = seq(rr[1], rr[2], length=tbreak)
    }
    nbreaks = length(tbreak)
    tbreak  = sort(tbreak)
    # We use the full tails if extended breaks are required
    if (ext) {
        tbreak[1] = -Inf
        tbreak[nbreaks] = Inf
    }
    y = table(cut(tstat, tbreak))

    # The function
    ff = function(param, df, nq, ng, y, tbreak) {
            delta = c(0, param[1:(nq-1)])
            beta = param[nq:(2*nq-2)]
            p = exp(beta)/(1+sum(exp(beta)))
            p = c(1-sum(p), p)
            FF = matrix(0, nrow=length(tbreak), ncol=nq)
            FF[,1] = pt(tbreak,df=df)
            for (j in 2:nq){
                FF[,j] = pt(tbreak,df=df, ncp=delta[j])
            }
            ff = apply(FF,2,diff)
            fit = ng* c(ff %*% p)
            ## -sum(log(fit)*y) ## modified to the following line by Long Qu
            -sum(log( pmax(fit, min(fit[fit>0])))*y)
    }

    # Prepare the starting values; constraints are handled by transformation
    param0 = c(delta, log(p1)-log(p0))
    oo = optim(param0, ff, df=df, nq=nq, ng=n, y=y, tbreak=tbreak, ...)
    
    # Nice output
    p1 = oo$par[nq:(2*nq-2 )]
    p1 = exp(p1)/(1+sum(exp(p1)))    
    p0 = 1-sum(p1)
    delta = oo$par[1:nq-1]
    AIC = 2*oo$value + 2*length(delta)
    # Only components with absolute delta greater than a given threshold 
    # contribute to the estimation proper!
    p0.est = p0 + sum(p1[abs(delta)<threshold.delta])
    #list(p0.est=p0.est, p0.raw=p0, p1=p1, D = delta/sf, delta = delta,          threshold.delta=threshold.delta, AIC=AIC, opt=oo) ## modified to the following line by Long Qu
    ans=list(p0.est=p0.est, p0.raw=p0, p1=p1, D = delta/sf, delta = delta,          threshold.delta=threshold.delta, AIC=AIC, opt=oo,	 					 data=list(tstat=tstat, df=n1+n2-2), pi0=p0)
    class(ans)='discTMix'
    ans
}



fitted.discTMix=#fitted.values.discTMix=
function(object, ...)
{
    if(any(is.na(object))) return (NA)
    dtmat=matrix(dt(object$data$tstat, object$data$df), length(object$data$tstat), length(object$p1)+1)
    for(i in 1:length(object$p1)+1)
        dtmat[,i]=suppressWarnings(dt(object$data$tstat, object$data$df, object$delta[i-1]))
    drop(dtmat%*%c(object$p0.raw, object$p1))
}
