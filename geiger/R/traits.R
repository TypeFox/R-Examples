fitContinuous=function(
phy,
dat,
SE = 0,
model=c("BM", "OU", "EB", "trend", "lambda", "kappa", "delta", "drift", "white"),
bounds=list(),
control=list(method=c("subplex","L-BFGS-B"), niter=100, FAIL=1e200, hessian=FALSE, CI=0.95),
ncores=NULL, ...
)
{

    # SE: can be vector or single value (numeric, NULL, or NA); vector can include NA
    # opt: a list with elements 'method', 'niter', 'FAIL'; 'method' may include several optimization methods
    # bounds: a list with elements specifying constraint(s): e.g., bounds=list(alpha=c(0,1))
    # control: a list with elements specifying method to compute likelihood
    # ...: node state dataframe

    # data matching
	td=treedata(phy, dat)
        ## add check to make sure only unique data used
        if (nrow(td$data) != length(unique(rownames(td$data))))
            stop("Multiple records per tip label")

	phy=td$phy
	dat=td$data
	dd=dim(dat)
	trts=dd[2]
	if(trts>1){
		nm=colnames(dat)
		res=lapply(1:trts, function(idx){
            fitContinuous(phy, dat[,idx], SE=SE, model=model, bounds=bounds, control=control)
        })
		names(res)=nm
        class(res)=c("gfits", class(res))
		return(res)
	} else {
		dat=dat[,1]
	}


    # CONTROL OBJECT for optimization
	ct=list(method=c("subplex","L-BFGS-B"), niter=100, FAIL=1e200, hessian=FALSE, CI=0.95)
	if(any(!names(control)%in%names(ct)->tmp)) warning("Unexpected 'control' parameters:\n\t", paste(names(control)[tmp], collapse="\n\t"), sep="")
	control=control[which(!tmp)]
	if("method"%in%names(control)) control$method=match.arg(control$method, c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent", "subplex"), several.ok=TRUE)
	ct[names(control)]=control
	if(ct$niter<2) stop("'niter' must be equal to or greater than 2")
	ct$hessian_P=1-ct$CI

	model=match.arg(model, c("BM", "OU", "EB", "trend", "lambda", "kappa", "delta", "drift", "white"))

    # CONTROL OBJECT for likelihood
        if(model=="OU" & !is.ultrametric(phy)){
            warning("Non-ultrametric tree with OU model, using VCV method.")
            #con=list(method="vcv",backend="R")
            #con[names(control)]=control
            lik=ou.lik(phy, dat, SE, model, ...)
        } else {
            con=list(method="pruning",backend="C")
            con[names(control)]=control
            lik=bm.lik(phy,dat,SE,model,...)
        }
        attr(lik, "model")=model
	argn=argn(lik)

    ## CONSTRUCT BOUNDS ##
	mn=c(-500, -500, (log(10^(-5))/max(branching.times(phy))), -100, -100, -500, -500, -500, -500)
	mx=c(100, 1, -0.000001, 100, 100, 0, 0, log(2.999999), 100)
	bnds=as.data.frame(cbind(mn, mx))
	bnds$typ=c("exp", "exp", "nat", "nat", "nat", "exp", "exp", "exp", "exp")
	rownames(bnds)=c("sigsq", "alpha", "a", "drift", "slope", "lambda", "kappa", "delta", "SE")
	bnds$model=c("BM", "OU", "EB", "drift", "trend", "lambda", "kappa", "delta", "SE")
	typs=bnds[argn, "typ"]

    # User bounds
	if(length(bounds)>0){
		mm=match(names(bounds), rownames(bnds))
		if(any(is.na(mm))){
			warning("Unexpected 'bounds' parameters:\n\t", paste(names(bounds)[is.na(mm)], collapse="\n\t"), sep="")
		}
		mm=mm[!is.na(mm)]

		if(length(mm)){
			for(i in 1:length(mm)){
				ww=mm[i]
				tmp=sort(bounds[[i]])
				if(bnds$typ[ww]=="exp") {
					if(any(tmp==0)) tmp[tmp==0]=exp(-500)
					bnd=log(tmp)
				} else {
					bnd=tmp
				}
				bnds[ww,c("mn","mx")]=bnd
			}
		}
	}
    if(any(!is.finite(as.matrix(bnds[,c("mn", "mx")])))) {
        stop("All bounds should be finite")
    }

	par=argn[1]


    ## likelihood function for optimizer (with modified space)
	xx=function(p){
		pars=ifelse(typs=="exp", exp(p), p)
		tmp=-lik(pars, root="max")
		if(is.infinite(tmp)) tmp=ct$FAIL
        if(is.na(tmp)) tmp=ct$FAIL
		tmp
	}

    # boxconstrain (from diversitree) to avoid out-of-bounds values
	boxconstrain=function (f, lower, upper, fail.value)
	{
		function(x) {
			if (any(x < lower | x > upper)) fail.value else f(x)
		}
	}
	f=boxconstrain(xx, bnds[argn,"mn"], bnds[argn,"mx"], fail.value=ct$FAIL)


    ## STARTING POINT -- problematic models ##
	if(par%in%c("alpha","lambda","delta","kappa")){
		bmstart=try(.bm.smartstart(phy,dat),silent=TRUE)
		if(inherits(bmstart, "try-error")) bmstart=0.01
        bmstart=log(bmstart)
	}
    #    rt=max(abs(dat))

    ## OPTIMIZATION ##
	mm=matrix(NA, nrow=ct$niter, ncol=length(argn)+2)
	mt=character(ct$niter)
    min=bnds[argn,"mn"]
    max=bnds[argn,"mx"]

    # mclapply or lapply
    fxopt=.get.parallel(ncores)

    # 'method' optimization
	out=fxopt(1:ct$niter, function(i){
		bnds$st=sapply(1:nrow(bnds), function(x) runif(1, bnds$mn[x], bnds$mx[x]))
		start=bnds[argn,"st"]

        ## OU ##
		if(par=="alpha"){
            oustart=log(.ou.smartstart(dat, unlist(exp(bnds["alpha",c("mx","mn")]))))
			if(i==1 | runif(1)<0.25) start[match(c("sigsq", par), argn)]=c(bmstart, oustart)
			if(runif(1) < 0.5) start[match(c("sigsq", par), argn)]=c(0,oustart)
		}

        ## PAGEL MODELS ##
		if(par%in%c("lambda","delta","kappa")){
			ww=match(par, rownames(bnds))
			if(runif(1)<0.5){
				if(runif(1)<0.5){
					start[match(c("sigsq", par), argn)]=c(bmstart, bnds$mx[ww])
				} else {
					start[match(c("sigsq", par), argn)]=c(bmstart, bnds$mn[ww])
				}
			}
		}

        ## WHITE NOISE ##
		if(par=="white"){
			if(runif(1)<0.5){
				start[match("sigsq", argn)]=var(dat)
			}
		}
		names(start)=argn

        # resolve method
		if(length(argn)==1) {
			method="Brent"
		} else {
            method=sample(ct$method,1)
		}

		if(method=="subplex"){
			op<-try(suppressWarnings(subplex(par=start, fn=f, control=list(reltol = .Machine$double.eps^0.25, parscale = rep(0.1, length(argn))), hessian=ct$hessian)),silent=TRUE)
		} else {
			op<-try(suppressWarnings(optim(par=start, fn=f, upper=max, lower=min, method=method, hessian=ct$hessian)),silent=TRUE)
		}

		if(!inherits(op,"try-error")){
            op$method=method
			op$value=-op$value
			names(op)[names(op)=="value"]="lnL"
			names(op$par)=argn
			op$par=sapply(1:length(typs), function(x) if(typs[x]=="exp") return(exp(op$par[x])) else return(op$par[x]))
			op$mm=c(op$par, op$lnL, op$convergence)
		} else {
			op=list(par=structure(rep(NA, length(argn)), names=argn), lnL=-Inf, convergence=1, method="FAIL")
		}

        op                                #return(op)
	})

    for(i in 1:length(out)){
        cur=out[[i]]
        if(cur$method!="FAIL") mm[i,]=cur$mm
        mt[i]=cur$method
    }

	res=mm
	colnames(res)=c(argn, "lnL","convergence")
	rownames(res)=mt

    ## HANDLE OPTIMIZER OUTPUT ##
	colnames(mm)=c(argn, "lnL", "convergence")
	conv=mm[,"convergence"]==0
	mm=mm[,-which(colnames(mm)=="convergence")]
	valid=apply(mm, 1, function(x) !any(is.na(x)))
	if(sum(valid & conv)>=1){
		mm=matrix(mm[valid,], nrow=sum(valid), dimnames=dimnames(mm))
		mt=mt[valid]
		out=out[valid]
		mm=mm[z<-min(which(mm[,"lnL"]==max(mm[,"lnL"]))),]
	} else {
		z=NA
		mm=c(rep(NA, length(argn)), -Inf)
		names(mm)=c(argn,"lnL")
	}
	zz=mm[-which(names(mm)%in%c("lnL"))]
	mm=as.list(mm)
    tmp=lik(unlist(mm[argn]))
    mm=c(mm[argn], z0=attributes(tmp)$ROOT.MAX, mm[names(mm)[!names(mm)%in%argn]])

	mm$method=ifelse(is.na(z), NA, mt[z])
	mm$k=length(argn)+1

    ## HESSIAN-based CI of par estimates
	if(ct$hessian){
		hessian=out[[z]]$hessian
		CI=.bnd.hessian(hessian, zz, typs, ct$hessian_P)
		if(!all(is.na(CI))){
			if(is.constrained(lik)){
				CI=rbind(lik(CI[1,], pars.only=TRUE), rbind(lik(CI[2,], pars.only=TRUE)))
			}
			dimnames(hessian)=NULL
			rownames(CI)=c("lb", "ub")
		}
	} else {
		hessian=NULL
		CI=NULL
	}


    # check estimates against bounds #
	range=as.data.frame(cbind(min, max))
	range$typ=typs
	range$mn=ifelse(range$typ=="exp", exp(range$min), range$min)
	range$mx=ifelse(range$typ=="exp", exp(range$max), range$max)
	par=mm[argn]
	rownames(range)=argn
	chk=sapply(1:length(par), function(idx){
        p=par[[idx]]
        if(!is.na(p)){
            return((p<=range$mn[idx] | p>=range$mx[idx]))
        } else {
            return(FALSE)
        }
    })
	if(any(chk)){
		warning(paste("Parameter estimates appear at bounds:\n\t", paste(names(par)[chk], collapse="\n\t", sep=""), sep=""))
	}

	mm=.aic(mm, n=length(dat))

    # RETURN OBJECT
	mm$CI=CI
	mm$hessian=hessian
    res=list(lik=lik, bnd=range[,c("mn", "mx")], res=res, opt=mm)
    class(res)=c("gfit", class(res))
    return(res)
}

.slot.attribute=function(x, attb, pos=1L){
	if(x%in%attb) {
		return(attb)
	} else {
		y=character(length(attb)+length(x))
		x.idx=c(pos:(pos+length(x)-1L))
		attb.idx=(1:length(y))[-x.idx]
		y[x.idx]=x
		y[attb.idx]=attb
		return(y)
	}
}



.fix.root.bm=function(root, cache){
    rtidx=cache$root
    cache$y["m",rtidx]=root
    attr(cache$y,"given")[rtidx]=as.integer(TRUE)
    cache
}

bm.lik=function(phy, dat, SE = NA, model=c("BM", "OU", "EB", "trend", "lambda", "kappa", "delta", "drift", "white"), ...){
	model=match.arg(model, c("BM", "OU", "EB", "trend", "lambda", "kappa", "delta", "drift", "white"))

    cache=.prepare.bm.univariate(phy, dat, SE=SE, ...)
    cache$ordering=attributes(cache$phy)$order ## SHOULD BE POSTORDER
    cache$N = cache$n.tip
    cache$n = cache$n.node
    cache$nn = (cache$root+1):(cache$N+cache$n)
    cache$intorder = as.integer(cache$order[-length(cache$order)])
    cache$tiporder = as.integer(1:cache$N)
    cache$z = length(cache$len)

	# function for reshaping edges by model
	FUN=switch(model,
        BM=.null.cache(cache),
        OU=.ou.cache(cache),
        EB=.eb.cache(cache),
        trend=.trend.cache(cache),
        lambda=.lambda.cache(cache),
        kappa=.kappa.cache(cache),
        delta=.delta.cache(cache),
        drift=.null.cache(cache),
        white=.white.cache(cache)
    )

    ll.bm.direct=function(cache, sigsq, q=NULL, drift=NULL, se=NULL){
        n.cache=cache

        given=attr(n.cache$y,"given")

        ## q
        if(is.null(q)) {
            llf=FUN()
        } else {
            llf=FUN(q)
        }
        ll=llf$len

        ## drift
        dd=0
        if(!is.null(drift)) dd=drift

        ## se
        adjvar=as.integer(attr(n.cache$y,"adjse"))
        adjSE=any(adjvar==1)
		.xxSE=function(cache){
            vv=cache$y["s",]^2
            ff=function(x){
				if(any(adjvar==1->ww)){
					vv[which(ww)]=x^2
					return(vv)
				} else {
					return(vv)
				}
			}
			return(ff)
		}
		modSE=.xxSE(n.cache)
		vv=as.numeric(modSE(se))

        ## PARAMETERS
        datC=list(
            len = as.numeric(ll),
            intorder = as.integer(n.cache$intorder),
            tiporder = as.integer(n.cache$tiporder),
            root = as.integer(n.cache$root),
            y = as.numeric(n.cache$y["m", ]),
            var = as.numeric(vv),
            n = as.integer(n.cache$z),
            given = as.integer(given),
            descRight = as.integer(n.cache$children[ ,1]),
            descLeft = as.integer(n.cache$children[, 2]),
            drift = as.numeric(dd)
        )
        #print(datC)
        parsC=as.numeric(rep(sigsq, n.cache$z))

        out = .Call("bm_direct", dat = datC, pars = parsC, package = "geiger")
        loglik <- sum(out$lq)
		if(is.na(loglik)) loglik=-Inf
        attr(loglik, "ROOT.MAX")=out$initM[datC$root]
        class(loglik)=c("glnL", class(loglik))
        return(loglik)
    }
    class(ll.bm.direct) <- c("bm", "dtlik", "function")

    ## EXPORT LIKELIHOOD FUNCTION
    fx_exporter=function(){

        ## OPTIONAL arguments
        attb=c()

        if(!is.null(qq<-argn(FUN))){
            adjQ=TRUE
            attb=c(attb, qq)
        } else {
            adjQ=FALSE
        }

        #sigsq
        attb=c(attb, "sigsq")

        #SE
        if(any(attr(cache$y, "adjse")==1)) {
            attb=c(attb, "SE")
        }

        #drift
        if(model=="drift") {
            attb=c(attb, "drift")
        }

        cache$attb=attb ## current attributes (potentially modified with 'recache' below)

        lik <- function(pars, ...) {

            ## ADJUSTMENTS of cache
            recache=function(nodes=NULL, root="max", cache){
                r.cache=cache
                if(root=="max"){
                    rtmx=TRUE
                } else if(root%in%c("obs", "given")){
                    rtmx=FALSE
                    r.cache$attb=c(cache$attb, "z0")
                } else {
                    stop("unusable 'root' type specified")
                }
                r.cache$ROOT.MAX=rtmx

                if(!is.null(nodes)) {
                    m=r.cache$y["m",]
                    s=r.cache$y["s",]
                    g=attr(r.cache$y, "given")
                    nn=r.cache$nn
                    r.cache$y=.cache.y.nodes(m, s, g, nn, r.cache$phy, nodes=nodes)
                }
                r.cache
            }
            rcache=recache(..., cache=cache)
            attb=rcache$attb

            if(missing(pars)) stop(paste("The following 'pars' are expected:\n\t", paste(attb, collapse="\n\t", sep=""), sep=""))

            pars=.repars(pars, attb)
            names(pars)=attb

            if(adjQ) q = pars[[qq]] else q = NULL
            sigsq = pars[["sigsq"]]
            if("SE"%in%attb) se=pars[["SE"]] else se=NULL
            if("drift"%in%attb) drift=-pars[["drift"]] else drift=0
            if("z0"%in%attb) rcache=.fix.root.bm(pars[["z0"]], rcache)

            ll = ll.bm.direct(cache=rcache, sigsq=sigsq, q=q, drift=drift, se=se)
            return(ll)
        }
        attr(lik, "argn") = attb
        attr(lik, "cache") <- cache
        class(lik) = c("bm", "function")
        lik
    }
    likfx=fx_exporter()
    return(likfx)
}

ou.lik <- function (phy, dat, SE = NA, ...)
{
    model="OU"
    if(is.na(SE)){SE=0}
    cache = .prepare.bm.univariate(phy, dat, SE = SE)
    cache$ordering = attributes(cache$phy)$order
    cache$N = cache$n.tip
    cache$n = cache$n.node
    cache$nn = (cache$root + 1):(cache$N + cache$n)
    cache$intorder = as.integer(cache$order[-length(cache$order)])
    cache$tiporder = as.integer(1:cache$N)
    cache$z = length(cache$len)
    ll.ou.vcv <- function(cache, alpha, sigsq, z0, se){
        N <- cache$N
        phyvcv <- vcv.phylo(cache$phy)
        ouV <- .ou.vcv(phyvcv)
        V <- sigsq*ouV(alpha)
        if(!is.null(se)){
            diag(V) <- diag(V)+se^2
        } else {
          if(any(attr(cache$y, "given")==1)){
            var <- cache$y['s',][attributes(cache$y)$given==1]^2
            var[is.na(var)] <- 0
            diag(V) <- diag(V)+var
          }
        }
        ci <- solve(V)
        invV <- ci
        if(is.null(z0)){
            o <- rep(1, length(cache$dat))
            m1 <- 1/(t(o) %*% ci %*% o)
            m2 <- t(o) %*% ci %*% cache$dat
            mu <- m1 %*% m2
        } else {mu <- z0}
        logd <- determinant(V, logarithm=TRUE)
        logl <- -0.5 * (t(cache$dat-mu) %*% invV %*%
                        (cache$dat-mu)) - 0.5 * logd$modulus -
                            0.5 * (N * log(2 * pi))
        attributes(logl) <- NULL
        if (is.na(logl)) logl = -Inf
        attr(logl, "ROOT.MAX") = mu
        class(logl) = c("glnL", class(logl))
        return(logl)
    }
    class(ll.ou.vcv) <- c("bm","dtlik","function")
    fx_exporter = function() {
        attb=c("alpha","sigsq")
            cache$attb <- attb
            if (any(attr(cache$y, "adjse") == 1)) {
                attb = c(attb, "SE")
            }
        lik <- function(pars, root="max", ...){
            attb=c("alpha","sigsq")
            cache$attb <- attb
            if (any(attr(cache$y, "adjse") == 1)) {
                attb = c(attb, "SE")
            }
            if (root == "max"){
                rtmx = TRUE
            } else {
                if (root %in% c("obs", "given")) {
                    attb = c(attb,"z0")
                } else {
                    stop("unusable 'root' type specified")
                }
            }
            if (missing(pars))
                stop(paste("The following 'pars' are expected:\n\t",
                           paste(attb, collapse = "\n\t", sep = ""), sep = ""))
            pars = .repars(pars, attb)
            names(pars) = attb
            if ("alpha" %in% attb)
                alpha = pars[["alpha"]]
            sigsq = pars[["sigsq"]]
            if ("SE" %in% attb){
                se = pars[["SE"]]
            } else se = NULL
            if ("z0" %in% attb) {
                z0 = pars[["z0"]]
            } else {z0 = NULL}
            ll = ll.ou.vcv(cache = cache, alpha=alpha, sigsq = sigsq,
                z0=z0, se = se)
            return(ll)
        }
        attr(lik, "argn") = attb
        attr(lik, "cache") <- cache
        class(lik) = c("bm", "function")
        lik
    }
    likfx <- fx_exporter()
    return(likfx)
}

.reshape.constraint.m=function(m){
	k=unique(dim(m))
	if(length(k)>1) stop("'m' must be a square matrix")
	diag(m)=NA
	tt=table(m)
	map=cbind(as.integer(names(tt)), seq_along(1:max(as.integer(length(tt)))))
	z=match(m, map[,1])
	m[]=map[z,2]
	class(m)=c("constraint.m",class(m))
	m
}


## MAIN FUNCTION for OPTIMIZATION
# to replace geiger:::fitContinuous

fitDiscrete=function(
phy,
dat,
model=c("ER","SYM","ARD","meristic"),
transform=c("none", "EB","lambda", "kappa", "delta", "white"),
bounds=list(),
control=list(method=c("subplex","L-BFGS-B"), niter=100, FAIL=1e200, hessian=FALSE, CI=0.95),
ncores=NULL,
...)
{

## NOTE: 'model' can be a constraint matrix
#
#		transform="none"
#		model="SYM"
#		bounds=list()
#		control=list(hessian=TRUE)

	td=treedata(phy, dat)

        ## add check to make sure only unique data used
        if (nrow(td$data) != length(unique(rownames(td$data))))
            stop("Multiple records per tip label")

	phy=td$phy
	dat=td$data
	dd=dim(dat)
	trts=dd[2]
	if(trts>1){
		nm=colnames(dat)
		res=lapply(1:trts, function(idx){
				   fitDiscrete(phy, dat[,idx], model=model, transform=transform, bounds=bounds, control=control, ...)
				   })
		names(res)=nm
        class(res)=c("gfits", class(res))
		return(res)
	} else {
		tmp=dat[,1]
		names(tmp)=rownames(dat)
		dat=tmp
		#if(!all(is.integer(dat))) stop("supply 'dat' as a vector (or matrix) of positive integers")

		# this changes the discrete data to 1:n and remembers the original charStates
		fdat<-as.factor(dat)
		charStates<-levels(fdat)
		k<-nlevels(fdat)

		ndat<-as.numeric(fdat)
    	names(ndat) <- names(fdat)

	}



	constrain=model
	model=match.arg(transform, c("none", "EB", "lambda", "kappa", "delta", "white"))

# CONTROL OBJECT
	ct=list(method=c("subplex","L-BFGS-B"), niter=ifelse(model=="none", 50, 100), FAIL=1e200, hessian=FALSE, CI=0.95)
	if(any(!names(control)%in%names(ct)->tmp)) warning("Unexpected 'control' parameters:\n\t", paste(names(control)[tmp], collapse="\n\t"), sep="")
	control=control[which(!tmp)]
	if("method"%in%names(control)) control$method=match.arg(control$method, c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent", "subplex"), several.ok=TRUE)
	ct[names(control)]=control
	if(ct$niter<2) stop("'niter' must be equal to or greater than 2")
	ct$hessian_P=1-ct$CI


	lik=mkn.lik(phy, ndat, constrain=constrain, transform=model, ...)
    attr(lik, "transform")=model
	if(model=="white") return(list(opt=lik))
	argn=unlist(argn(lik))

	# translation of original character states
	attr(lik, "levels")<-charStates


## CONSTRUCT BOUNDS ##
#	minTrans<-log(1)-log(sum(phy$edge.length))
    minTrans<--500
	maxTrans<- log(100*Nedge(phy))-log(sum(phy$edge.length))
	mn=c(-10, -500, -500, -5, minTrans)
	mx=c(10, 0, 0, log(2.999999), maxTrans)
	bnds=as.data.frame(cbind(mn, mx))
	bnds$typ=c("nat", "exp", "exp", "exp", "exp")
	rownames(bnds)=c("a", "lambda", "kappa", "delta", "trns")
	bnds$model=c("EB", "lambda", "kappa", "delta", "trns")

	parnm=ifelse(argn%in%rownames(bnds), argn, "trns")
	typs=bnds[parnm, "typ"]

# User bounds
	if(length(bounds)>0){
		mm=match(names(bounds), rownames(bnds))
		if(any(is.na(mm))){
			warning("Unexpected 'bounds' parameters:\n\t", paste(names(bounds)[is.na(mm)], collapse="\n\t"), sep="")
		}
		mm=mm[!is.na(mm)]
		if(length(mm)){
			for(i in 1:length(mm)){
				ww=mm[i]
				tmp=sort(bounds[[i]])
				if(bnds$typ[ww]=="exp") {
					if(any(tmp==0)) tmp[tmp==0]=exp(-500)
					bnd=log(tmp)
				} else {
					bnd=tmp
				}
				bnds[ww,c("mn","mx")]=bnd
			}
		}
	}
    if(any(!is.finite(as.matrix(bnds[,c("mn", "mx")])))) {
        stop("All bounds should be finite")
    }


	par=argn[1]

## likelihood function for optimizer (with modified space)
	xx=function(p){
		pars=ifelse(typs=="exp", exp(p), p)
		tmp=-lik(pars, root="obs", root.p=NULL)
		if(is.infinite(tmp)) tmp=ct$FAIL
        if(is.na(tmp)) tmp=ct$FAIL
		tmp
	}



# boxconstrain from diversitree
	boxconstrain=function (f, lower, upper, fail.value)
	{
		function(x) {
			if (any(x < lower | x > upper)) fail.value else f(x)
		}
	}
	f=boxconstrain(xx, min, max, fail.value=ct$FAIL)

# transition rates
	smarttrns=log(seq(max(c(exp(bnds["trns","mn"]),0.01)), min(c(exp(bnds["trns","mx"]),2)), length.out=10))

## OPTIMIZATION ##
	mm=matrix(NA, nrow=ct$niter, ncol=length(argn)+2)
	mt=character(ct$niter)
    min=bnds[parnm,"mn"]
    max=bnds[parnm,"mx"]

    # mclapply or lapply
    fxopt=.get.parallel(ncores)

    # 'method' optimization
	out=fxopt(1:ct$niter, function(i){
		bnds$st=sapply(1:nrow(bnds), function(x) runif(1, bnds$mn[x], bnds$mx[x]))
		start=bnds[parnm,"st"]

        ## PAGEL MODELS ##
		if(par%in%c("lambda","delta","kappa")){
			ww=match(par, rownames(bnds))
			if(runif(1)<0.5){
				if(runif(1)<0.5){
					start[match(par, argn)]=bnds$mx[ww]
				} else {
					start[match(par, argn)]=bnds$mn[ww]
				}
			}
		}

        ## TRANSITION RATES
		if(runif(1)<0.5){
			start[which(parnm=="trns")][]=smarttrns[sample(1:length(smarttrns),1)]
		}
		names(start)=argn

        # resolve method
		if(length(argn)==1) {
			method="Brent"
		} else {
            method=sample(ct$method,1)
		}

		if(method=="subplex"){
			op<-try(suppressWarnings(subplex(par=start, fn=f, control=list(reltol = .Machine$double.eps^.25, parscale = rep(0.1, length(argn))), hessian=ct$hessian)),silent=TRUE)
		} else {
			op<-try(suppressWarnings(optim(par=start, fn=f, upper=max, lower=min, method=method, hessian=ct$hessian)),silent=TRUE)
		}

		if(!inherits(op,"try-error")){
            op$method=method
			op$value=-op$value
			names(op)[names(op)=="value"]="lnL"
			names(op$par)=argn

			op$par=sapply(1:length(typs), function(x) if(typs[x]=="exp") return(exp(op$par[x])) else return(op$par[x]))
			op$mm=c(op$par, op$lnL, op$convergence)
		} else {
			op=list(par=structure(rep(NA, length(argn)), names=argn), lnL=-Inf, convergence=1, method="FAIL")
		}
        return(op)
    })

    for(i in 1:length(out)){
        cur=out[[i]]
        if(cur$method!="FAIL") mm[i,]=cur$mm
        mt[i]=cur$method
    }

	res=mm
	colnames(res)=c(argn,"lnL","convergence")
	rownames(res)=mt

## HANDLE OPTIMIZER OUTPUT ##
	colnames(mm)=c(argn, "lnL", "convergence")
	conv=mm[,"convergence"]==0
	mm=mm[,-which(colnames(mm)=="convergence")]
	valid=apply(mm, 1, function(x) !any(is.na(x)))
	if(sum(valid & conv)>=1){
		mm=matrix(mm[valid,], nrow=sum(valid), dimnames=dimnames(mm))
		mt=mt[valid]
		out=out[valid]
		mm=mm[z<-min(which(mm[,"lnL"]==max(mm[,"lnL"]))),]
		mod=mt[z]
	} else {
		mod=NA
		z=NA
		mm=c(rep(NA, length(argn)), -Inf)
		names(mm)=c(argn,"lnL")
	}
	k=length(argn)
	llx=which(names(mm)=="lnL")
	ll=mm[[llx]]
	mm=mm[-llx]

## HESSIAN-based CI of par estimates
	if(ct$hessian){
#		print(mm)
#		print(partp)
#		print(out[[z]]$hessian)
		hessian=out[[z]]$hessian
		CI=.bnd.hessian(hessian, mm, typs, ct$hessian_P)
		if(!all(is.na(CI))){
			if(is.constrained(lik)){
				CI=rbind(lik(CI[1,], pars.only=TRUE), rbind(lik(CI[2,], pars.only=TRUE)))
			}
			dimnames(hessian)=NULL
			rownames(CI)=c("lb", "ub")
		}
	} else {
		hessian=NULL
		CI=NULL
	}

# resolve all transition estimates (if constrained model)
	trn=!names(mm)%in%c("a", "lambda", "kappa", "delta")
	if(is.constrained(lik)){
		constr=TRUE
		allq=lik(mm, pars.only=TRUE)
	} else {
		constr=FALSE
		allq=mm[argn(lik)]
	}
	allq=allq[!names(allq)%in%names(mm)[!trn]]

	qparnm=rep("trns", length(allq))
	allparnm=c(parnm[!trn],qparnm)
	min=bnds[allparnm,"mn"]
	max=bnds[allparnm,"mx"]
	typs=bnds[allparnm, "typ"]
	argn=c(argn[!trn], names(allq))
	mmcomplete=c(mm[!trn], allq, ll)
	names(mmcomplete)=c(argn, "lnL")

	mm=as.list(mmcomplete)
	mm$method=mod
	mm$k=k
	mm=.aic(mm, n=length(dat))


# check estimates against bounds #
	range=as.data.frame(cbind(min, max))
	range$typ=typs
	range$mn=ifelse(range$typ=="exp", exp(range$min), range$min)
	range$mx=ifelse(range$typ=="exp", exp(range$max), range$max)
	par=mm[argn]
	rownames(range)=argn
	chk=sapply(1:length(par), function(idx){
			   p=par[[idx]]
			   if(!is.na(p)){
			   return((p<=range$mn[idx] | p>=range$mx[idx]))
			   } else {
			   return(FALSE)
			   }
			   })
	if(any(chk)){
		warning(paste("Parameter estimates appear at bounds:\n\t", paste(names(par)[chk], collapse="\n\t", sep=""), sep=""))
	}

# RETURN OBJECT
	mm$CI=CI
	mm$hessian=hessian

	res=list(lik=lik, bnd=range[,c("mn", "mx")], res=res, opt=mm)
    class(res)=c("gfit", class(res))
    return(res)
}

## WORKHORSE -- built from diversitree:::make.mkn
# likelihood function creation
mkn.lik=function(
phy,
dat,
constrain=c("ER","SYM","ARD","meristic"),
transform=c("none", "EB", "lambda", "kappa", "delta", "white"),
...)
{
    phy=reorder(phy, "postorder")

    # control object for make.mkn()
	ct = list(method="exp")
	if(ct$method!="exp") stop(paste("method",sQuote(ct$method),"is not currently supported",sep=" "))

    # primary cache
	k<-nlevels(as.factor(dat))
    if(is.character(dat)) dat=structure(as.factor(dat), names=names(dat))
    if(is.factor(dat)){
        levels=levels(dat)
        dat=structure(as.integer(dat), names=names(dat))
    } else {
        levels=sort(unique(dat))
    }
    if(k==2) if(all(constrain=="SYM")) constrain="ER"
    control <- .check.control.mkn(ct, k)
    cache <- .make.cache.mkn(phy, dat, k, strict=TRUE, control=ct, ultrametric=FALSE)
	cache$ordering=attributes(cache$info$phy)$order

    # tree transforms
	trns=match.arg(transform, c("none", "EB","lambda", "kappa", "delta", "white"))

	FUN=switch(trns,
    none=.null.cache(cache),
    EB=.eb.cache(cache),
    lambda=.lambda.cache(cache),
    kappa=.kappa.cache(cache),
    delta=.delta.cache(cache),
    white=white.mkn(cache$states))

	if(trns=="white") return(FUN)

	ll.mkn=function(cache, control, ...) {
		k <- cache$info$k
		f.pars <- .make.pars.mkn(k)
		f.pij <- .make.pij.mkn(cache$info, control)
		idx.tip <- cache$idx.tip
		n.tip <- cache$n.tip
		n <- length(cache$len)
		map <- t(sapply(1:k, function(i) (1:k) + (i - 1) * k))
		idx.tip <- cbind(c(map[cache$states, ]), rep(seq_len(n.tip), k))
		children.C <- .toC.int(t(cache$children))
		order.C <- .toC.int(cache$order)

		.ll.mkn.exp=function(q, pars, intermediates=FALSE, preset = NULL) { # based on diversitree:::make.all.branches.mkn.exp
			if(is.null(argn(FUN))) new=FUN() else new=FUN(q)

			len.uniq <- sort(unique(new$len))
			len.idx <- match(new$len, len.uniq)

			if (!is.null(preset)) stop("Preset values not allowed")
			pij <- f.pij(len.uniq, pars)[, len.idx]
			lq <- numeric(n)
			branch.init <- branch.base <- matrix(NA, k, n)
			storage.mode(branch.init) <- "numeric"
			ans <- matrix(pij[idx.tip], n.tip, k)
			q <- rowSums(ans)
			branch.base[, seq_len(n.tip)] <- t.default(ans/q)
			lq[seq_len(n.tip)] <- log(q)
			ans <- .C("r_mkn_core", k = as.integer(k), n = length(order.C) -
            1L, order = order.C, children = children.C, pij = pij,
            init = branch.init, base = branch.base, lq = lq,
            NAOK = TRUE, PACKAGE="geiger")


			list(init = ans$init, base = ans$base, lq = ans$lq, vals = ans$init[, cache$root], pij = pij)
		}

        # build likelihood function
		attb=c(argn(FUN), cache$info$argn)
        rt=function(root="obs", root.p=NULL){
                    return(list(root=root, root.p=root.p))
        }
		if(is.null(argn(FUN))){ # NO TRANSFORM
			ll=function(pars, ...){
                rx=rt(...)
				qmat=f.pars(pars)
				ans=.ll.mkn.exp(q=NULL, pars=qmat, intermediates=FALSE)
				.rootfunc.mkn(ans, qmat, root=rx$root, root.p=rx$root.p, intermediates=FALSE)
			}
		} else {
			ll=function(pars, ...){ # TREE TRANSFORM
                rx=rt(...)
				qmat=f.pars(pars[-1])
				ans=.ll.mkn.exp(q=pars[1], pars=qmat, intermediates=FALSE)
				.rootfunc.mkn(ans, qmat, root=rx$root, root.p=rx$root.p, intermediates=FALSE)
			}

		}
		class(ll) <- c("mkn", "dtlik", "function")
		attr(ll,"argn") <- attb
        if(!is.null(levels)) attr(ll, "levels")=levels
		return(ll)
	}

	tmp=ll.mkn(cache, control)

    ## CONSTRAINTS
	if(!all(constrain=="ARD")){
		if(is.character(constrain)){
			cc=match.arg(constrain, c("ER","SYM","ARD","meristic"))
			tmp=.constrain.k(tmp, model=cc, ...)
		} else {
			if(is.matrix(constrain)){
				if(ncol(constrain)==max(dat)){
					tmp=.constrain.m(tmp, m=constrain)
				}
			} else {
				stop("'constrain' must be supplied as a dummy matrix representing constraints on transition classes")
			}
		}
	}
	lik=function(pars, ...){
        if(missing(pars)) stop(paste("The following 'pars' are expected:\n\t", paste(argn(tmp), collapse="\n\t", sep=""), sep=""))

		pars=.repars(pars, argn(tmp))
		tmp(pars, ...)
	}
	attributes(lik)<-attributes(tmp)
    attr(lik, "trns")<-!argn(lik)%in%argn(FUN)
	lik
}

aov.phylo=function(formula, phy, nsim=1000, test=c("Wilks", "Pillai", "Hotelling-Lawley", "Roy"), ...){
    xx=lapply(all.vars(formula), get)
    flag="'formula' must be of the form 'dat~group', where 'group' is a named factor vector and 'dat' is a data matrix or named vector"

    if(!is.factor(xx[[2]])) stop(flag)
    if(is.null(names(xx[[2]]))) stop(flag)

    yy=merge(xx[[1]], xx[[2]], by=0)
    if(nrow(yy)==0) stop(flag)
    rownames(yy)=yy[,1]
    yy=yy[,-1]

    tmp<-treedata(phy, yy, sort=TRUE)
    phy=tmp$phy
    yy=yy[phy$tip.label,]

    group=structure(yy[,ncol(yy)], names=rownames(yy))
    dat=as.matrix(yy[,-ncol(yy)])
    rownames(dat)=rownames(yy)

    s<-ratematrix(phy, dat)

    multivar=ifelse(ncol(dat)>1, TRUE, FALSE)
    if(multivar){
        test=match.arg(test, c("Wilks", "Pillai", "Hotelling-Lawley", "Roy"))
        m=summary.manova(mod<-manova(dat~group), test=test)
        f.data=m[[4]][1,2]
        FUN=function(xx) summary.manova(manova(as.matrix(xx)~group), test=test)[[4]][1,2]
        sims<-sim.char(phy, s, nsim=nsim)
        f.null<-apply(sims, 3, FUN)
        out=as.data.frame(m[[4]])
        attr(out, "heading")=c("Multivariate Analysis of Variance Table\n","Response: dat")
    } else {
        test=NULL
        m=anova(mod<-lm(dat~group))
        f.data<-m[1,4]
        FUN=function(xx) anova(lm(xx~group))[1,4]
        out=as.data.frame(m)
    }
    colnames(out)=gsub(" ", "-", colnames(out))
    sims<-sim.char(phy, s, nsim=nsim)
    f.null<-apply(sims, 3, FUN)
    if(multivar) {
    	if(test=="Wilks") {
    		p.phylo = (sum(f.null < f.data) + 1)/(nsim + 1)
    	} else {
    		p.phylo = (sum(f.null > f.data) + 1)/(nsim + 1)
    	}
    } else {
    	p.phylo = (sum(f.null > f.data) + 1)/(nsim + 1)
    }
    out$'Pr(>F) given phy'=c(p.phylo, NA)
    class(out) <- c("anova", "data.frame")
    print(out, ...)
    attr(mod, "summary")=out
    return(mod)
}
