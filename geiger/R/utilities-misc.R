.geigerwarn <- function(...) warning("the called function is currently in development and is not fully vetted", ...);

coef.gfit <- function(object, ...) {
    if (is.constrained(object$lik)) p=names(object$lik(argn(object$lik),pars.only=TRUE)) else p=argn(object$lik)
    if ("bm"%in%class(object$lik)) p=c(p, "z0");
    unlist(object$opt[p]);
}

coef.gfits <- function(object, ...) {
    lapply(object, coef);
}

logLik.gfit <- function(object, ...) {
    object$opt$lnL;
}

logLik.gfits <- function(object, ...) {
    lapply(object, function(x) x$opt$lnL);
}

# get aic-weights
# x is a named vector of AIC scores
aicw <- function (x) {
    if (!inherits(x, "numeric"))
        stop("aic scores must be of class 'numeric'")

	aic <- x;
	best <- min(aic);
	delta <- aic - best;
	sumDelta <- sum(exp(-0.5 * delta));
	w <- (exp(-0.5 * delta)/sumDelta);

	results <- data.frame(fit = aic, delta = delta, w = w);
	rownames(results) <- names(aic);

	results
}

#general printing utility for ensuring equal numbers of characters within columns and defining spacing between columns
#author: JM EASTMAN 2010
#note: works only for numeric dataframes
.print.table=function(df,digits=4,buffer=5){
	if(length(buffer) != ncol(df) | length(buffer)==1) buffer=rep(buffer[1],ncol(df))
	if(length(digits) != ncol(df) | length(digits)==1) digits=rep(digits[1],ncol(df))
	ss=sapply(round(df),nchar)
	lar=df>1
	nn=sapply(names(df),nchar)

    # find longest string
	strw=sapply(1:ncol(df), function(x) max(nn, max(1,(ss[lar])+digits[x],na.rm=TRUE),na.rm=TRUE))
	pr.df=data.frame(do.call(cbind, lapply(1:ncol(df), function(x) sprintf(paste("%",(strw[x]+buffer[x]),".",digits[x],"f",sep=""),df[,x]))))
	names(pr.df)=names(df)
	rownames(pr.df)=rownames(df)
	print(pr.df)
}

# ooh, this is nice.
.get.parallel <- function (ncores = NULL, ...)
{
    if ((Sys.getenv("R_PARALLEL") == "FALSE")) {
        fx <- function(X, FUN, ...) lapply(X, FUN, ...)
    }
    else {
        if (.check.parallel() & Sys.info()["sysname"] != "Windows") {
            if (is.null(ncores)) {
                ncores <- parallel::detectCores()
            }
            fx <- function(X, FUN, ...) parallel::mclapply(X, FUN, ...,
                mc.silent = TRUE, mc.cores = ncores)
        }
        else {
            fx <- function(X, FUN, ...) lapply(X, FUN, ...)
        }
    }
    return(fx)
}

.check.parallel <- function() {
	if (.gui.check()) {
		return (FALSE);
	}
	tmp <- rownames(installed.packages());
	if ("parallel" %in% tmp) {
		return(TRUE);
	} else {
		return(FALSE);
	}
}

# prevent parallel from loading if gui
.gui.check <- function () {
	if (!is.na(Sys.getenv()["R_GUI_APP_VERSION"])) {
		return (TRUE);
	} else {
		return (FALSE);
	}
}

.transparency <- function (col, alpha) {
    tmp <- col2rgb(col)/255;
    rgb(tmp[1, ], tmp[2, ], tmp[3, ], alpha = alpha);
}

.withinrange <- function (x, min, max) {
    a <- sign(x - min);
    b <- sign(x - max);
    if (abs(a + b) == 2) {
		return(FALSE);
    } else {
    	return(TRUE);
    }
}

.basename.noext <- function(path="") {
	return(sub("[.][^.]*$", "", basename(path), perl=TRUE));
}

# make this more general. joe shmo won't expect v to contain lik + k; split'em up
# .aic <- function (lnL, n, k) {
	# res <- NULL;
	# res$aic <- 2 * k - 2 * lnL;
	# res$aicc <- 2 * k * (n/(n - k - 1)) - (2 * lnL);
	# return(res);
#}

# v: has object 'lnL' and 'k'
# .aic <- function (v, n) {
	# res <- NULL;
	# res$aic <- 2 * v$k - 2 * lnL;
	# res$aicc <- 2 * v$k * (n/(n - v$k - 1)) - (2 * v$lnL);
	# return(res);
# }



.aic <- function (v, n) {
# v: has object 'lnL' and 'k'
	v$aic <- 2 * v$k - 2 * v$lnL;
    #v$aicc <- 2 * v$k * (n - 1)/(n - v$k - 2) - 2 * v$lnL; # wrong
    v$aicc <- 2 * v$k * (n/(n - v$k - 1)) - (2 * v$lnL);
	return (v);
}

.resolve.executable=function(package="geiger"){
	packagedir=system.file(package=package)
	execs=lapply(d<-dir(paste(packagedir,"exec",sep="/")), function(x) {paste(packagedir, "exec", x, sep="/")})
	names(execs)=.basename.noext(d)
	return(execs)
}


#rjmcmc utility for initiating a proposal width for Markov sampling
#author: JM EASTMAN 2010
#modified: 02.26.2011 to use spline() in computation of best prop.width
#deprecates calibrate.proposalwidth
calibrate.rjmcmc <- function(phy, dat, nstep=10000, widths=2^(-3:3), type=c("bm",
"rbm", "jump-bm", "jump-rbm"), ...) {
	model=match.arg(type, c("bm",
    "rbm", "jump-bm", "jump-rbm"))

    acceptance.rates=sapply(widths, function(x) rjmcmc.bm(phy=phy, dat=dat, ngen=nstep, samp=1, prop.width=x, summary=FALSE, type="rbm", ...)$acceptrate)

	aa=sapply(acceptance.rates, .withinrange, 0.20, 0.80)
	df=data.frame(width=widths, acceptrate=acceptance.rates)
	.print.table(df)

	if(any(aa)){
		acc=acceptance.rates[aa]
		wid=widths[aa]
		s=spline(log(widths,base=2),acceptance.rates)
		choice=2^(mean(s$x[s$y==max(s$y)]))
	} else {
		stop("No proposal width found with acceptance rates between 0.20 and 0.80.")
	}

	return(choice)
}

.fix.rjmcmc.matrix=function(mat, colnames){
	if(any(duplicated(colnames(mat))) | any(duplicated(colnames)) | any(is.na(colnames)) | any(is.na(colnames(mat)))){
		stop("Non-unique hash keys encountered as node indentifiers")
	}
	mm=matrix(NA, nrow=nrow(mat), ncol=length(colnames))
	vv=match(colnames(mat)->nm,colnames)
	new=as.matrix(mat[,!is.na(vv)])
	colnames(new)=colnames(mat)[!is.na(vv)]
	mm[,match(colnames(new), colnames)]=new
	attrb=attributes(mat)
	attrb=attrb[-which(names(attrb)%in%c("dim","dimnames"))]
	attributes(mm)[names(attrb)]=attrb
	colnames(mm)=colnames
	mm
}


#load=function(x, ...){
#	tmp=all(sapply(x, function(y) file.info(y)$isdir==TRUE))
#	type=ifelse(tmp, "dir", "rda")
#	if(type=="dir"){
#		return(load.rjmcmc(x, ...))
#	} else {
#		dots=list(...)
#		if("envir"%in%names(dots)){
#			return(load(x, envir=dots$envir))
#		} else {
#			return(load(x, envir=parent.frame()))
#		}
#	}
#}

#merges samples from multiple independent Markov chains (generated by rjmcmc.bm() or mcmc.levy())
#author: JM EASTMAN 2010
load.rjmcmc <- function(x, phy=NULL, burnin = NULL, thin = NULL, ...){
    #	single x & no tree [currently returning original tree]
    #	single x & tree [currently returning original tree with hashes based on 'phy' and 'hashtips']
    #	many x & tree [currently returning samples based on 'phy' and 'hashtips', summarized on 'phy']

	dirs=x

	z=list(...)
	if("hashtips"%in%names(z)) hashtips=z$hashtips else hashtips=NULL

	if(length(dirs)==1) {
		return(.subset.auteurRAW(get(load(paste(x, dir(x, pattern="samples.rda"), sep="/"))), burnin=burnin, thin=thin, phy=phy, hashtips=hashtips))
	}
	FUN=.get.parallel()
	raw<-FUN(dirs, function(x) get(load(paste(x, dir(x, pattern="samples.rda"), sep="/"))))

	# TREES: collect trees and resolve variable topologies
	trees=lapply(raw, "[[", "phy")
	class(trees)="multiPhylo"
	uu=unique(trees)
	if(length(uu)>1){
		if(is.null(phy)) stop("Encountered multiple topologies: a summary tree is necessary (supplied via the 'phy' argument)")
	} else {
		if(!is.null(phy)) {
			tmp=unique(list(uu,phy))
			if(length(tmp)==1) {
				phy=phy
			}
		} else {
			phy=uu[[1]]
		}
	}
	class(phy)="phylo"
	if(is.null(hashtips)) hashtips=phy$tip.label
	phy=hashes.phylo(phy, hashtips)

	# SAMPLES: collect individual runs: construct matrices for all unique edges (only necessary if multiple trees run, differing in topology)
	samples<-FUN(raw, function(x) {
        y=.subset.auteurRAW(x, burnin=burnin, thin=thin, phy=phy, hashtips=hashtips)
        for(i in 1:length(y)){
            if(!is.null(attributes(y[[i]])$parm)) y[[i]]=.fix.rjmcmc.matrix(y[[i]], phy$hash[-c(Ntip(phy)+1)])
        }
        y$tree=x$phy

        y})

	# TREES
	trees=lapply(samples, "[[", "tree")
	class(trees)="multiPhylo"

	for(i in 1:length(samples)){
		samples[[i]]$tree=NULL
	}

	# CODA: convert matrices to mcmc.list
	csamp=suppressWarnings(to.coda(samples))
	csamp$trees=trees

	res=to.auteur(csamp, phy=phy)
	res
}

.detachDiversitree=function(){
	tt=try(detach(package:diversitree), silent=TRUE)
	if(!inherits(tt, "try-error")) warning("'diversitree' functions have been masked; use 'require(diversitree)' to reload the package")

}


## USE INSTEAD OF samples$edger (replacing .edgewise.rjmcmc)
## CREATE post-hoc when given a tree (if given a tree)
.edger=function(samples, phy=NULL, hashtips=NULL){
	orig=samples$phy
	class(orig)="phylo"
	if(is.null(phy)) phy=orig
	class(phy)="phylo"
	if(is.null(hashtips)) hashtips=phy$tip.label
	orig=hashes.phylo(orig, hashtips)

	desc=.cache.descendants(orig)
	rootd=desc$fdesc[[Ntip(orig)+1]]

	cache=list(desc=desc,phy=orig)

	nd=orig$edge[,2]
	vv=numeric(length(nd))

	prep.string=function(string, sep="\t"){
		st=as.numeric(unlist(strsplit(string, split=sep)))
		nds=st[seq(1,length(st),by=2)]
		val=st[seq(2,length(st), by=2)]

		return(list(nodes=nds, values=val))
	}

	to.vector=function(nodes, values, heritable=FALSE){
		d=duplicated(nodes)
		nodes=nodes[!d]
		values=values[!d]
		if(heritable){
			for(i in 1:length(nodes)) vv=.assigndescendants(vv, nodes[i], values[i], exclude=nodes, cache)
		} else {
			vv[match(nodes, nd)]=values
		}
		vv
	}

	vd=function(string, sep="\t", par=c("shifts", "jumps"), tree.only=FALSE){

		if(tree.only) return(orig)
		mm<-ss<-matrix(0, nrow=length(string), ncol=length(nd))

		par=match.arg(par, c("jumps", "shifts"))
		if(par%in%c("jumps")) {
			heritable=FALSE
		} else if(par%in%c("shifts")) {
			heritable=TRUE
		} else {
			stop("'par' not recognized")
		}

		for(k in 1:length(string)){
			cur=prep.string(string[k], sep="\t")
			mm[k,]=to.vector(cur$nodes, cur$values, heritable=heritable)
			if(heritable==TRUE){ # log shift points for heritable parameter
				cur$nodes=cur$nodes[-c(1:length(rootd))]
				if(length(cur$nodes)){
					cur$values=rep(1,length(cur$nodes))
					ss[k,]=to.vector(cur$nodes, cur$values, heritable=FALSE)
				}
			}
		}
		colnames(mm)=orig$hash[nd]
		attr(mm, "hashtips")=hashtips
		if(heritable){
			attr(mm, "shifts")=ss
		}
		mm
	}
	vd
}

.subset.auteurRAW=function(x, burnin=NULL, thin=NULL, phy=NULL, hashtips=NULL){

    ## returns rjmcmc output based on 'phy' and (or) 'hashtips'

	.detachDiversitree()

	trim=function(obj, burnin=NULL, thin=NULL){
		orig="matrix"
		if(!"matrix"%in%class(obj)){
			orig=class(obj)
			obj=as.matrix(obj)
		}
		nn=nrow(obj)
		if(!is.null(burnin)){
			if(!.withinrange(burnin, 0, 1)) stop("Supply 'burnin' as a fraction (between 0 and 1)")
			if(burnin>0){
				bb=ceiling(burnin*nn)
				obj=as.matrix(obj[-c(1:bb),])
			}
		}
		if(!is.null(thin)){
			if(!.is.wholenumber(thin)) stop("Supply 'thin' as an integer.")
			if(thin>1){
				tt=seq(1, nrow(obj), by=thin)
				obj=as.matrix(obj[tt,])
			}

		}
		if(orig!="matrix") obj=as(obj,orig)
		return(obj)
	}

	samples=x

    ## LOG FILE ##
	logf=read.table(samples$log, header=TRUE)
	st=trim(storig<-logf$state, burnin, thin)
	tmp=diff(st)
	logf=logf[match(st, storig),]
	logf$ppos=logf[,"lnL"]+logf[,"lnLp"]
	logf=logf[,!colnames(logf)%in%c("lnLp", "qlnL.p", "qlnL.h")]
	if(length(thinned<-unique(tmp))!=1) stop("Encountered uninterpretable log file.")
	trace=.coda(data=logf[,-which(names(logf)=="state")], start=min(st), end=max(st), thin=thinned)
	class(trace)=c("rjmcmc", class(trace), class(unclass(trace)))
	mcpar<-mcparx<-attr(trace,"mcpar")
	names(mcparx)=c("start","end","thin")


    ## EDGER
	if(!is.null(phy)) {
		class(phy)="phylo"
		if(is.null(hashtips)) hashtips=phy$tip.label
	} else {
		phy=samples$phy
	}
	if(!is.null(hashtips)) hashtips=hashtips else hashtips=samples$phy$tip.label

    #	phy=hashes.phylo(phy, hashtips)

	FUN=.edger(samples, phy, hashtips)

    ## PARMS
	util=c("hasher", "edger", "phy", "log")
	pars=names(samples)[!names(samples)%in%util]
	mats=lapply(pars, function(x) {
        cur=trim(samples[[x]], burnin, thin)
        z=FUN(cur, sep="\t", par=x)
        #				colnames(z)=attr(z, "hash")[attr(z,"nodes")]
        #				attr(z, "hash")=NULL
        thin=mcparx[["thin"]]
        ngen=mcparx[["end"]]
        attr(z,"parm")=TRUE
        ff=mcmc(data=z, start=thin, thin=thin)
        attr(ff, "mcpar")=mcpar
        class(ff)=c("rjmcmc", class(ff), class(unclass(ff)))
        #				ff=.fix.rjmcmc.matrix(ff, phy$hash[-c(Ntip(phy)+1)])
        ff
    })
	names(mats)=pars
	mats$rates=mats$shifts
	sft=attr(mats$shifts, "shifts")
	attributes(sft)=attributes(mats$rates)
	mats$shifts=sft
	attr(mats$shifts,"shifts")<-attr(mats$rates,"shifts")<-NULL

	mats$log=trace
	mats$phy=FUN(tree.only=TRUE)
	class(mats)=c("auteurMCMC", class(mats))

    #	if(!is.null(phy)) mats$phy=hashes.rjmcmc(mats, phy) else mats$phy=hashes.rjmcmc(mats, samples$phy)

	return(mats)
}

.coda=function (data = NA, start = 1, end = numeric(0), thin = 1) ## from coda:::mcmc
{
    if (is.matrix(data)) {
        niter <- nrow(data)
        nvar <- ncol(data)
    }
    else if (is.data.frame(data)) {
 	    if (!all(dd<-apply(data, 2, function(x) all(is.na(x))))) {
 	    	data=data[,which(!dd)]
            #          stop("Data frame contains non-numeric values")
		}
        data <- as.matrix(data)
        niter <- nrow(data)
        nvar <- ncol(data)
    }
    else {
        niter <- length(data)
        nvar <- 1
    }
    thin <- round(thin)
    if (length(start) > 1)
	stop("Invalid start")
    if (length(end) > 1)
	stop("Invalid end")
    if (length(thin) != 1)
	stop("Invalid thin")
    if (missing(end))
	end <- start + (niter - 1) * thin
    else if (missing(start))
	start <- end - (niter - 1) * thin
    nobs <- floor((end - start)/thin + 1)
    if (niter < nobs)
	stop("Start, end and thin incompatible with data")
    else {
        end <- start + thin * (nobs - 1)
        if (nobs < niter)
		data <- data[1:nobs, , drop = FALSE]
    }
    attr(data, "mcpar") <- c(start, end, thin)
    attr(data, "class") <- "mcmc"
    data
}



## CONVERT OUTPUT to CODA object (mcmc.list)
to.coda=function(obj){

	to.mcmc.list=function(obj){
        # expect auteurMCMCMC
		nrun=attr(obj,"nrun")
		mcpar=attr(obj, "mcpar")
		ll=mcpar[2]
		rr=unique(rownames(obj))
		if(length(rr)!=nrun) stop("Encountered corrupted 'obj'.")
		samples=list()
		for(i in 1:length(rr)){
			ww=which(rownames(obj)==rr[i])
			tmp=obj[ww,]
			rownames(tmp)=NULL
			aa=attr(obj, "mcpar")
			aa[2]=ll
			attr(tmp,"mcpar")=aa
			class(tmp)=c("rjmcmc","mcmc",class(tmp))
			samples[[i]]=tmp
		}
		res=mcmc.list(samples)
		class(res)=c("rjmcmc.list", class(res))
		return(res)
	}

	if("auteurMCMCMC"%in%class(obj)){
		nn=names(obj)
		nn=nn[!nn%in%c("phy","trees")]
		phy=obj$phy
		trees=obj$trees
		res=lapply(obj[nn], to.mcmc.list)
		names(res)=nn
		class(res)=c("codaMCMCMC", class(res))
		res$phy=phy
		res$trees=trees
		return(res)
	}

	if("auteurMCMC"%in%class(obj)) {
		return(obj)
	}

	if("rjmcmcmc"%in%class(obj)){
		return(to.mcmc.list(obj))
	}

	if("rjmcmc"%in%class(obj)) {
		return(obj)
	}

    # ELSE: expect list of 'auteurMCMC' objects
	zz=sapply(obj, function(x) "auteurMCMC"%in%class(x))
	if(!all(zz)) stop("Supply 'obj' as a list of 'auteurMCMC' objects -- see to.auteur()")
	warning("'load.rjmcmc()' is recommended for combining multiple runs.")
	samples=obj
	nn=unique(unlist(lapply(samples, names)))

    # ensure identity of trees
	trees=sapply(samples, function(x) x$phy)
	class(trees)="multiPhylo"
	uu=unique(trees)
	if(length(uu)==1){
		phy=uu[[1]]
	} else {
		phy=NULL
	}
	nn=nn[nn!="phy"]

	samples=list()
	for(i in nn) {
		tmp=lapply(1:length(obj), function(x) obj[[x]][[i]])
		cur=mcmc.list(tmp)
		samples[[which(nn==i)]]=cur
	}
	names(samples)=nn
	samples$phy=phy
	class(samples)=c("codaMCMCMC", class(samples))
	samples
}


## CONVERT OUTPUT to rjMCMC object
to.auteur=function(obj, phy=NULL, ...){
	samples=obj
	if("auteurRAW"%in%class(samples)) return(.subset.auteurRAW(samples, phy=phy, ...))

	zz=list(...)
	if(any(c("burnin","thin")%in%names(zz))){
		warning("'burnin' and 'thin' have no effect in this context: use load.rjmcmc() on the original output directories.")
	}

	if("mcmc.list"%in%class(samples)){
		res=.intercalate.rjmcmc(samples)
		if(!is.null(phy)) warning("returning intercalated samples ('phy' has no effect in this context)")
		return(res)
	}

	if("codaMCMCMC"%in%class(samples)){
		nn=names(samples)
		excl=c("phy","trees", "hasher")
		nn=nn[!nn%in%excl]
		res=lapply(samples[nn], .intercalate.rjmcmc)
		names(res)=nn
		class(res)=c("auteurMCMCMC", class(res))
		excl=excl[excl%in%names(samples)]
		if(length(excl)) res[excl]=samples[excl]
		if(!is.null(phy)) res$phy=hashes.rjmcmc(res, phy)
		return(res)
	}

	if(any(c("auteurMCMC","auteurMCMCMC")%in%class(samples))){
		if(!is.null(phy)){
			samples$phy=hashes.rjmcmc(samples, phy)
			return(samples)
		} else {
			return(samples)
		}
	}

	if(any(c("rjmcmcmc","rjmcmc")%in%class(samples))){
		if(!is.null(phy)) warning("returning unmodified samples ('phy' has no effect in this context)")
		return(samples)
	}
}


#utility for intercalating a list of matrices into a single 'intercalated' sample.
#expects 'mcmc.list' as 'obj'
#author: JM EASTMAN 2010
.intercalate.rjmcmc <-
function(obj) {
	list.obj=obj
	if(!"mcmc.list"%in%class(obj)) stop("Supply 'obj' as an 'mcmc.list' object.")
	if(length(list.obj)<2) warning("Nothing to be done... too few objects in list.")

	vv=sapply(list.obj, is.vector)
	if(all(vv)) list.obj=lapply(list.obj, function(x) {y=matrix(x,nrow=1); colnames(y)=names(x); return(y)})

	exclude.attr=c("shifts", "nodes")
	zz=sapply(list.obj, function(x) {
        y=names(attributes(x)->aa)
        digest(aa[!y%in%exclude.attr])
    })

    # resolve attributes for rjmcmcmc object
	if(length(unique(zz))==1){
		aa=attributes(list.obj[[1]])
		dim=aa$dim
		dimnames=aa$dimnames
		attrb=aa[!names(aa)%in%c("dim","dimnames","class",exclude.attr)]
	} else {
		stop("Chains in 'obj' do near appear comparable.")
	}

    # construct intercalated matrix
	n=length(list.obj)
	indices=lapply(1:n, function(x) seq(x, dim[1]*n, by=n))

	out.array=array(dim=c(dim[1]*n, dim[2]))
	for(nn in 1:n){
		out.array[indices[[nn]],]=as.matrix(list.obj[[nn]])
	}

	rownames(out.array)=rep(1:n, length=nrow(out.array))
	colnames(out.array)=dimnames[[2]]

	attributes(out.array)[names(attrb)]=attrb
	attr(out.array, "nrun")=n
	class(out.array)=c("rjmcmcmc", "mcmc", class(unclass(out.array)))
	return(out.array)
}


hashes.rjmcmc=function(obj, phy){
	checktree=function(hashes, hashtips, phy){
		if(!"hphylo"%in%class(phy)) phy=hashes.phylo(phy, hashtips)
		phyh=phy$hash[phy$edge[,2]]
		if(any(phyh%in%hashes)) return(phy) else stop("'phy' cannot be matched to 'obj'.")
        #		if(all(phyh%in%hashes)) return(phy) else stop("'phy' cannot be matched to 'obj'.")
	}

	if(any(c("rjmcmcmc","rjmcmc")%in%class(obj))){ # dealing with single parameter (or log)
		if(!is.null(attr(obj, "parm"))) {
			ht=attr(obj,"hashtips")
			hh=colnames(obj)
			return(checktree(hh, ht, phy))
		} else {
			stop("'phy' has no effect in this context")
		}
	}

	if(!any(c("auteurMCMCMC","auteurMCMC")%in%class(obj))) stop("Supply an object of class 'auteurMCMCMC' or 'auteurMCMC'")
	nn=names(obj)[sapply(obj, function(x) if(!is.null(attr(x,"parm"))) return(TRUE) else return(FALSE))]
	par=obj[nn]

	ht=c()
	hh=c()
	for(i in 1:length(par)){
		cur=attr(par[[i]],"hashtips")
		curh=colnames(par[[i]])
		if(!length(ht)) {
			ht=cur
			hh=curh
		} else {
			if(!all(ht==cur)){
				stop("Parameters in 'obj' have unexpected colnames.")
			}
			if(!all(hh%in%curh)){
				stop("Parameters in 'obj' have unexpected colnames.")
			}
		}
	}
	return(checktree(hh, ht, phy))
}

## RUN CONTROL ##
make.gbm=function(phy, dat, SE=NA, type=c("bm", "rbm", "jump-bm", "jump-rbm"), ...){
	type=match.arg(type, c("rbm","bm","jump-rbm","jump-bm"))
	if(type%in%c("bm", "rbm", "jump-bm", "jump-rbm")){
		flavor="gbm"
	}
	switch(flavor,
    gbm=.prepare.gbm.univariate(phy, dat, SE=SE, type=type, ...)
    )
}

.prepare.gbm.univariate=function(phy, dat, SE=NA, type=c("bm", "rbm", "jump-bm", "jump-rbm"), ...){

	type=match.arg(type, c("rbm","bm","jump-rbm","jump-bm"))

	con=list(
    method="direct",												# likelihood method: direct - pruning algorithm
    rate.lim=list(min=0, max=Inf),									# limits on rates
    root.lim=list(min=-10^3, max=10^3),							# limits on root
    se.lim=list(min=0, max=Inf),									# limits on SE
    constrainSHIFT=FALSE,											# limit number of local rates (under *relaxedBM)
    constrainJUMP=FALSE,											# limit number of jumps (under jump*)
    dlnSHIFT=NULL,													# shift prior (function; arg 'x'; returns ln(p[x]))
    dlnJUMP=NULL,													# jump prior  (function; arg 'x'; returns ln(p[x]))
    dlnROOT=NULL,													# root prior (function; arg 'x'; returns ln(p[x]))
    dlnRATE=function(x) dexp(x, rate=1/(10^3), log=TRUE),			# rate prior  (function; arg 'x'; returns ln(p[x]))
    dlnPULS=function(x) dexp(x, rate=1/(10^3), log=TRUE),			# pulse prior  (function; arg 'x'; returns ln(p[x]))
    dlnSE=function(x) dexp(x, rate=1/(10^3), log=TRUE),			# se prior	(function; arg 'x'; returns ln(p[x]))
    jump.lim=1,													# limit on number of jumps per branch
    excludeSHIFT=c(),												# edges to exclude for shifts
    excludeJUMP=c(),												# edges to exclude for jumps
    bm.jump=0.5,													# proposal density of 'bm' updates versus 'jump' updates
    mergesplit.shift=0.5,											# proposal density for 'mergesplit' (merging or splitting of rate classes) versus proposals that 'shift' a rate class in tree
    tune.scale=0.65,												# proposal density for 'tune' a single rate class or 'scale' all rate classes
    slide.mult=0.25,												# proposal density for 'slide' (sliding window) proposals versus 'mult' (multiplier) proposals
    prob.dimension=0.65,											# proposals for dimensionality: involves 'bm.jump', 'mergesplit.shift'
    prob.effect=0.30,												# proposals for effects of process: 'slide.mult', 'tune.scale'
    prob.root=0.02,												# proposals for effects of process: 'slide.mult', 'tune.scale'
    prob.SE=0.03,													# proposals for effects of process: 'slide.mult', 'tune.scale'
    prop.width=1,													# proposal width for sliding window and multiplier proposals
    sample.priors=FALSE,											# sample from priors only
    simple.start=TRUE,												# start with single rate class and no jumps -- overruled if is.numeric(constrainK) and is.numeric(constrainJ)
    summary=TRUE,													# used with calibrate.rjmcmc
    beta=1,														# exponent (of likelihood) used to create power-posteriors (Xie et al. 2010, Syst Biol)
    filebase="result",												# output directory
    primary.parameter="shifts"										# INTERNAL: for logging purposes
    )

	if(missing(phy) & missing(dat)) {
		# print direction to user if controller() called without data
		cat(paste(toupper("\n*** default priors and control parameters are as follows ***"),"\n\n",sep=""))
		cat("\n\nNOTE: settings within the control object can be adjusted through the additional arguments (...) of make.gbm()\n\n")

		print(con)

		cat("\n\nNOTE: settings within the control object can be adjusted through the additional arguments (...) of make.gbm()\n\n")
		stop("'phy' and 'dat' must minimally be supplied to this function")
	}

	## USER control-object
	xusr=list(...)
	usr=xusr[names(xusr)%in%names(con)]
	con[names(usr)]=usr
	missed=names(xusr)[!names(xusr)%in%names(usr)]
	if(length(missed)) cat(paste(paste(sQuote(missed),collapse=" "), "not expected in control object", sep=" "))

	## FINALIZE control-object
	# create cache
	con$method=match.arg(con$method, c("direct"))
	lltmp=make.bm.relaxed(phy, dat, SE, method=con$method)
	cache=attributes(lltmp)$cache
	attr(lltmp,"cache")=NULL
	cache$edge.length.cs=cumsum(cache$edge.length)
	con$lik=lltmp

	# check (or create) prior distributions
	nn=length(cache$edge.length)
	if(is.null(con$dlnSHIFT)){
		mm=nn-length(con$excludeSHIFT)
		dv=dcount(0:(max(c(1,mm))-1), FUN=dtpois, min=0, max=max(c(1,mm))-1, lambda=log(2))
		if(is.numeric(con$constrainSHIFT)){
			dv=.dunifn(con$constrainSHIFT)
		}
		con$dlnSHIFT=dv
	}
	if(is.null(con$dlnJUMP)){
		mm=nn-length(con$excludeJUMP)
		dv=dcount(0:mm, FUN=dlunif, min=0, max=mm, dzero=0.5)
		if(is.numeric(con$constrainJUMP)){
			dv=.dunifn(con$constrainJUMP)
		}
		con$dlnJUMP=dv

	}
	if(is.null(con$dlnROOT)){
		con$dlnROOT=function(x) dunif(x, min=con$root.lim$min, max=con$root.lim$max, log=TRUE)	# root prior  (function; arg 'x'; returns ln(p[x]))
	}

	if(any(is.infinite(sapply(con$root.lim, con$dlnROOT)))) stop("Infinite prior probability returned at bounds 'root.lim'")

	con$dlnSHIFT=.check.prior(con$dlnSHIFT, count=TRUE)
	con$dlnJUMP=.check.prior(con$dlnJUMP, count=TRUE)
	con$dlnRATE=.check.prior(con$dlnRATE, count=FALSE)
	con$dlnROOT=.check.prior(con$dlnROOT, count=FALSE)
	con$dlnPULS=.check.prior(con$dlnPULS, count=FALSE)
	con$dlnSE=.check.prior(con$dlnSE, count=FALSE)

	# check exclude vectors
	if(length(con$excludeSHIFT)){
		if(!all(con$excludeSHIFT%in%cache$nodes)) stop("Some edges in 'excludeSHIFT' not found in 'phy'.")
	}
	if(length(con$excludeJUMP)){
		if(!all(con$excludeJUMP%in%cache$nodes)) stop("Some edges in 'excludeJUMP' not found in 'phy'.")
	}

	# resolve model flavor
	tmp=.check.gbm(type, con$constrainJUMP, con$constrainSHIFT, cache)
	con$model=tmp$model
	con$constrainJUMP=tmp$constrainJUMP
	con$constrainSHIFT=tmp$constrainSHIFT
	con$algo=tmp$algo

	# check constraints on model dimensionality
	if(is.numeric(con$constrainSHIFT)){
		if(!.withinrange(con$constrainSHIFT, 0, length(cache$nodes)-1)) stop("'constrainSHIFT' must lie between 0 and nrow(phy$edge)-1.")
		con$algo="mcmc"
	} else {
		con$algo="rjmcmc"
		con$constrainSHIFT=NULL
	}
	if(is.numeric(con$constrainJUMP)){
		if(!.withinrange(con$constrainJUMP, 0, length(cache$nodes))) stop("'constrainJUMP' must lie between 0 and nrow(phy$edge).")
	} else {
		con$constrainJUMP=NULL
	}

	# proposal distributions
	if(sum(attributes(cache$y)$adjse)==0) {
        con$prob.SE=0
    } else {
        if(con$prob.SE==0) stop("'prob.SE' should be non-zero")
    }
	proposals=c(dim=con$prob.dimension, effect=con$prob.effect, root=con$prob.root, se=con$prob.SE)
	if(con$method=="reml") proposals=proposals[c("dim","effect")]
	prop.cs=cumsum(proposals*(1/sum(proposals)))
	names.subprop<-c("mergesplit","rootstate","ratetune","moveshift","ratescale","movejump","incjump","decjump","jumpvar","SE")
	n.subprop<-n.subaccept<-rep(0,length(names.subprop))
	names(n.subprop)<-names(n.subaccept)<-names.subprop
	con$prop.cs=prop.cs
	con$n.subprop=n.subprop
	con$n.subaccept=n.subaccept

	# file handling
	if(con$summary){
		parmbase=paste(con$model, con$filebase, sep=".")
		if(!file.exists(parmbase)) dir.create(parmbase)
		errorlog=paste(parmbase,paste(con$model, con$filebase, con$algo, "errors.log",sep="."),sep="/")
		runlog=file(paste(parmbase,paste(con$filebase, con$algo, "log",sep="."),sep="/"),open='w+')
		parms=list(principal=list(shifts=NULL, jumps=NULL), gen=NULL, lnL=NULL, lnLp=NULL, qlnL.p=NULL, qlnL.h=NULL, jumpvar=NULL, SE=NULL, root=NULL)
		con$runlog=runlog
		con$errorlog=errorlog
		con$parmbase=parmbase
		parl=c("shifts","jumps")
		parlgs=lapply(parl, function(f) file(paste(parmbase,paste(f, "txt",sep="."),sep="/"),open='w+'))
		names(parlgs)=parl
		con$parlogs=parlgs
		.parlog.rjmcmc(init=TRUE, end=FALSE, parameters=parms, con)
		con$parlogger=.ests.parlog.rjmcmc(cache)
	}

	# starting point
	start=.startingpt.bm(con, cache)

	return(list(control=con, cache=cache, start=start))
}


.check.gbm=function(type=c("rbm","bm","jump-rbm","jump-bm"), constrainJUMP=FALSE, constrainSHIFT=FALSE, cache){

	# resolve model
	type=match.arg(type, c("bm","rbm","jump-rbm","jump-bm"))
	if(type%in%c("rbm", "bm")){
		if(!is.numeric(constrainJUMP)) constrainJUMP=0
	}
	if(type%in%c("jump-bm", "bm")){
		if(!is.numeric(constrainSHIFT)) constrainSHIFT=0
	}

    if(type%in%c("jump-bm", "jump-rbm")) .geigerwarn(immediate.=TRUE)

	mod=expand.grid(constrainSHIFT=c(TRUE,FALSE),constrainJUMP=c(TRUE,FALSE))
	rownames(mod)=c("jump-rbm","jump-bm","rbm","bm")

	check.mod=function(parm, null=0){
		if(is.numeric(parm)) return(parm!=null) else return(!parm)
	}
	tmp=c(shift=check.mod(constrainSHIFT, null=0), jump=check.mod(constrainJUMP, null=0))
	model=rownames(mod)[modzz<-apply(mod, 1, function(x) all(x==tmp))]

	full=c("jump-relaxedBM", "jump-BM", "relaxedBM", "BM")

	if(model!=type){
		stop("Arguments 'type', 'constrainSHIFT', and (or) 'constrainJUMP' are inconsistent: see 'cache()'")
	} else {
		usrmodel=full[which(rownames(mod)==model)]
	}

	# resolve 'constrainJUMP' 'constrainSHIFT' and 'algo'
	if(is.numeric(constrainSHIFT)){
		if(!.withinrange(constrainSHIFT, 0, length(cache$nodes)-1)) stop("'constrainSHIFT' must lie between 0 and nrow(phy$edge)-1.")
		algo="mcmc"
	} else {
		algo="rjmcmc"
		constrainSHIFT=NULL
	}
	if(is.numeric(constrainJUMP)){
		if(!.withinrange(constrainJUMP, 0, length(cache$nodes))) stop("'constrainJUMP' must lie between 0 and nrow(phy$edge).")
	} else {
		constrainJUMP=NULL
	}

	return(list(model=usrmodel, constrainJUMP=constrainJUMP, constrainSHIFT=constrainSHIFT, algo=algo))


}


.ests.parlog.rjmcmc=function(cache){
	nm=cache$phy$edge[,2]
	root=cache$root
	rootd=cache$desc$fdesc[[root]]
	rr=match(rootd, nm)
	estlog=function(delta, values=NULL){
		xx=which(delta!=0)
		nd=c(nm[rr], nm[xx])
		if(length(delta)!=length(nm)) stop(paste("Expecting 'delta' of length ", length(nm), ".", sep=""))
		if(is.null(values)){
			values=delta
		} else {
			if(!length(values)==length(delta)) stop("'delta' and 'values' are mismatched.")
		}
		vv=values[c(rr, xx)]
		ests=paste(nd, vv, collapse="\t", sep="\t")
		ests
	}
	estlog
}

.read.ests.rjmcmc=function(control, cache){
	parlogs=control$parlogs
	mod=control$model
	log=control$runlog
	base=control$parmbase

	res=list()
	if(length(parlogs)){
		for(i in 1:length(parlogs)){
			res[[i]]=readLines(parlogs[[i]])
		}
		res
	}
}


#general phylogenetic utility for rapid computation of the maximum-likelihood estimate of the Brownian motion rate parameter (unless there are polytomies, in which case the function wraps geiger:::fitContinuous)
#author: L REVELL 2010, LJ HARMON 2009, and JM EASTMAN 2010
.bm.rate.mle <-
function(phy, dat){
	phy=multi2di(phy, random=TRUE)
	n=Ntip(phy)
	ic=try(pic(dat,phy),silent=TRUE)
	if(!inherits(ic, "try-error")) {
		r=mean(ic^2)*((n-1)/n)
		return(r)
	} else {
		return(var(dat)/mean(phy$edge.length[phy$edge.length>0]))
	}
}

## SIMULATION UTILITY ##
.startingpt.bm <- function(control, cache) {

	.rates.simulation=function(phy, exclude=NULL){
		drp=phy$edge[,2]%in%exclude
		nm=phy$edge[!drp,2]
		nt=Ntip(phy)
		sim=function(N) {
			n=1
			if(N==1) return(NULL)
			shifts=c()
			while(n<N){
				shifts=c(shifts, nm[z<-sample(1:length(nm), 1)])
				nm=nm[-z]
				n=n+1
			}
			return(c(sort(shifts[shifts>nt]), shifts[shifts<=nt]))

		}
		return(sim)
	}

	phy=cache$phy
	dat=cache$dat

	rate=.bm.rate.mle(phy,dat)
	if(!.withinrange(rate, control$rate.lim$min, control$rate.lim$max)) rate=runif(1, control$rate.lim$min, control$rate.lim$max)
	root=mean(dat)
	if(!.withinrange(root, control$root.lim$min, control$root.lim$max)) {
		root=runif(1, control$root.lim$min, control$root.lim$max)
	} else {
		root=.proposal.slidingwindow(root, control$prop.width, control$root.lim)$v
	}
	if(sum(attributes(cache$y)$adjse)==0) {
		se=NA
	} else {
		se=runif(1)
		if(!.withinrange(se, control$se.lim$min, control$se.lim$max)) se=runif(1, control$se.lim$min, control$se.lim$max)
	}

	nd=argn(control$lik)$rates
	tmp=numeric(length(nd))

	# rates
	if(control$simple.start & !is.numeric(control$constrainSHIFT)){
		rts=1
	} else {
		rts=ifelse(is.null(control$constrainSHIFT), .rcount(1, control$dlnSHIFT), control$constrainSHIFT+1)
	}
	rs=.rates.simulation(phy, control$excludeSHIFT)
	rat=rs(rts)
	rr<-dd<-tmp
	rr[]=rate
	if(length(rat)){
		for(i in 1:length(rat)){
			cur=rat[i]
			curd=c(cur, cache$desc$ades[[cur]])
			md=match(cur, nd)
			dd[md]=1
			curv=.proposal.multiplier(rate, control$prop.width, control$rate.lim)$v
			mm=match(curd, nd)
			rr[mm]=curv
		}
	}

    # jumps
	if(control$simple.start & !is.numeric(control$constrainJUMP)){
		jps=0
	} else {
		jps=ifelse(is.null(control$constrainJUMP), .rcount(1, control$dlnJUMP), control$constrainJUMP)
	}
	nds=phy$edge[,2]
	nds=nds[!nds%in%control$excludeJUMP]
	if(!is.null(control$constrainJUMP)) {
		if(length(nds)<control$constrainJUMP) stop(paste("Too few branches are available for jumps with restrictions given by 'constrainJUMP' and 'excludeJUMP':\n\tonly", length(nds), "branches are unrestricted", sep=" "))
	}

	j=0
	jj=rep(0,nrow(cache$phy$edge))
	if(jps>0){
		nd=rep(nds, control$jump.lim)
		jmp=nds[sample(1:length(nds), size=jps, replace=FALSE)]
		if(length(jmp)){
			je=.jumps.edgewise(phy)
			jj=je(jmp)
		} else {
			stop("Encountered unexpected error attempting to initialize vector of jumps")
		}
	}

	jv=rate*10

    rttmp=.proposal.multiplier(mean(rr), control$prop.width, control$rate.lim)$v
    rtrt=.link.root(rootd=cache$desc$fdesc[[cache$n.tip+1]], rttmp, nd, dd, rr)

	return(list(root=root, rates=rr, delta=dd, rootrate=rtrt, jumps=jj, jumpvar=jv, se=se))
}


#logging utility used by rjmcmc
#author: JM EASTMAN 2010
.parlog.rjmcmc <- function(init=FALSE, end=FALSE, parameters, control) {
	if(control$summary){
		parms=parameters
		parmbase=control$parmbase
		runlog=control$runlog
		method=control$method
		parlogs=control$parlogs
		primpar=control$primary.parameter
		reml=method=="reml"

		# add eol to files if terminating run
		if(end==TRUE) {
			if(length(parlogs)) sapply(1:length(parlogs), function(x) {write("\n", file=parlogs[[x]], append=TRUE); close(parlogs[[x]])})
			close(control$runlog)
			return()
		}

		if(reml) {
			parms=parms[names(parms)!="root"]
			pT=0
		} else {
			pT=control$dlnROOT(as.numeric(parms[["root"]]))
		}
		parnames=names(parms)
		ii=match(c("gen"),parnames)
		general=match(c("lnL","lnLp","qlnL.p","qlnL.h"),parnames)
		target=match("principal",parnames)
		logpars=names(parlogs)
		princpars=parms$principal
		accessory=match(parnames[-c(ii,target,general)],parnames)

		if(init) {
			write.table(paste("state", "min", "max", "median", paste(names(princpars),collapse="\t"), paste(parnames[accessory],collapse="\t"), paste(parnames[general],collapse="\t"), sep="\t"), file=runlog, quote=FALSE, col.names=FALSE, row.names=FALSE)
		} else {
			pp=princpars[[primpar]][[primpar]]
			range.p=range(pp)
			median.p=median(pp)
            #			pR=sum(control$dlnRATE(pp))
			ss=sapply(princpars, function(x) sum(x[["delta"]]))
			names(ss)=names(princpars)

			# to log file
			if(ss[["jumps"]]==0) {
				parms$jumpvar=0
                #				pV=0
			} else {
                #				pV=control$dlnRATE(parms$jumpvar)
			}
            #			pJ=control$dlnJUMP(ss[["jumps"]])

            #			pS=control$dlnSHIFT(ss[["shifts"]])

            # compute prior
            #			parms[["lnLp"]]=sum(c(pT, pR, pV, pJ, pS))

			msg<-paste(parms[[ii]], sprintf("%.3f", range.p[1]), sprintf("%.3f", range.p[2]), sprintf("%.3f", median.p), paste(ss, collapse="\t"), paste(sprintf("%.3f", parms[accessory]), collapse="\t"), paste(sprintf("%.2f", parms[general]),collapse="\t"),sep="\t")
			write(msg, file=runlog, append=TRUE)

            # to parameter files
			if(length(logpars)){
				for(i in 1:length(logpars)){
					cur=princpars[[logpars[i]]]
					curlog=parlogs[[i]]
					curpar=names(parlogs)[i]
					if(curpar%in%names(cur)) values=cur[[curpar]] else values=NULL
					tmp=control$parlogger(delta=cur$delta, values=values)
					write(tmp, file=curlog, append=TRUE)
				}
			}
		}
	}
}


#rjmcmc run diagnosis, generating tallies of proposed and accepted updates by class of proposal mechanism
#author: JM EASTMAN 2010
.load.rjmcmc <- function(control) {
	n.accept=control$n.subaccept
	n.props=control$n.subprop
	prop.names=names(control$n.subprop)
	df=data.frame(cbind(proposed=n.props, accepted=n.accept, adoptrate=n.accept/n.props))
	rownames(df)=prop.names
	if(control$summary){
		cat("\n\n",rep(" ",10),toupper(" sampling summary"),"\n")

		if(any(is.na(df))) df[is.na(df)]=0
		.print.table(df, digits=c(0,0,4), buffer=6)

		cat("\n\n")
		control$runlog=NULL
		control$errorlog=NULL
		control$parlogs=NULL
		return(control)
	} else {
		return(list(proposals=df, acceptrate=sum(df$accepted)/sum(df$proposed)))
	}
}


#logging utility used by rjmcmc
#author: JM EASTMAN 2010
## FIXME
.error.rjmcmc <-
function(i, proposal, mod.cur, mod.new, lnR, lnPrior, lnHastings, errorLog) {
	if(!file.exists(file=errorLog)) {
		write(paste("gen", "proposal", "cur.lnL", "new.lnL", "lnR", "ln.p", "ln.h", sep="\t"), file=errorLog)
	}
	write(paste(i, sprintf("%s", proposal), sprintf("%.3f", mod.cur), sprintf("%.3f", mod.new), sprintf("%.3f", lnR), sprintf("%.3f", lnPrior), sprintf("%.3f", lnHastings), sep="\t"),  file=errorLog, append=TRUE)

}

#utility for converting text files to .rda to compress output from rjmcmc
#author: JM EASTMAN 2011
.cleanup.rjmcmc <- function(control, cache){
	if(control$summary){
		parlogs=control$parlogs
		parlocs=sapply(parlogs, function(x) summary(x)$description)
		logloc=summary(control$runlog)$description


		.parlog.rjmcmc(init=FALSE, end=TRUE, parameters=NULL, control=control)

		samples=lapply(parlocs, function(x) {
            y=readLines(x)
            y=y[y!=""]
        })
		names(samples)=names(parlogs)
        #		cache$phy=hashes.phylo(cache$phy, cache$hashtips)
        #		samples$edger=.edgewise.rjmcmc(control, cache)
        #		samples$hasher=function(phy) hashes.phylo(phy, cache$hashtips)
        #		prior=list(rate=control$dlnRATE, shift=control$dlnSHIFT, jump=control$dlnJUMP, root=control$dlnROOT)
        #		samples$prior=prior
		samples$phy=cache$phy
		samples$log=logloc
		class(samples)=c("auteurRAW", class(samples))

		save(samples, file=paste(control$parmbase,paste(control$filebase,"samples","rda",sep="."),sep="/"))
		for(i in 1:length(parlocs)) unlink(parlocs[i])
	}
	.load.rjmcmc(control)
}
