ess=function(x){
	log=x
	if(!"mcmc"%in%class(log)) stop("Supply 'log' as an 'mcmc' object")
	rm=c("lnL.p", "lnL.h")
	log=log[,which(!colnames(log)%in%rm)]
	apply(log, 2, effectiveSize)
}

sem=function(x){
	val=x[!is.na(x)]
    if(!length(val)) return(NA)
	sd(val)/sqrt(length(val))
}



#general statistical function for computing a highest density region (whose proportion is given by 'hpd')
#author: R FITZ-JOHN 2011 (from diversitree)

hdr <-
function (z, hpd = 0.95, lim = NULL)
{
	hdr.uniroot <-
	function (z, p = 0.95, lim = NULL) 
	{
		xx <- sort(c(lim, seq(min(z), max(z), length = 1024)))
		ez <- ecdf(z)
		f <- suppressWarnings(approxfun(ez(xx), xx))
		fit <- suppressWarnings(optimize(function(x) f(x + p) - f(x), c(0, 1 - p)))
		if (inherits(fit, "try-error") || is.na(fit$objective)) 
		stop("HDR interval failure")
		ci <- fit$min
		f(c(ci, ci + p))
	}
	
    ci <- try(hdr.uniroot(z, hpd, lim), silent = TRUE)
    if (inherits(ci, "try-error")) {
        warning("HDR falling back on quantile-based intervals")
        ci <- as.numeric(quantile(z, c((1 - hpd)/2, 1/2 + hpd/2)))
    }
    ci
}




## HARMONIC MEAN ESTIMATOR of marginal likelihood
hme=function(x, scale=c("lnL", "logL", "p")){
	scale=match.arg(scale, c("p", "logL", "lnL"))
	const=0
	y=x
	if(scale=="lnL") {
		const=max(x)
		y=exp(x-const)
	}
	if(scale=="logL") {
		const=max(x)
		y=10^(x-const) 
	}
	zz=length(y)/sum(1/y)
	if(scale=="lnL") zz=log(zz)
	if(scale=="logL") zz=log(zz, 10)
	zz=zz+const
	zz
}

.bf.interpreter=function(x){
	set=unique(c(as.character(x[,"m1"]), as.character(x[,"m2"])))
	
	foo=function(interp){
		p=interp
		z=ifelse(p=="bare mention", "~", ifelse(p%in%c("positive","substantial"), ">", ifelse(p%in%c("strong"), ">>", ">>>")))
		return(z)	
	}
	
	idx=numeric(length(set)-1)
	for(i in 1:(length(set)-1)){
		y1=set[i]
		y2=set[i+1]
		
		idx[i]=which(x[,"m1"]==y1 & x[,"m2"]==y2)
	}
	res=sapply(x[idx,"interp"], foo)
	txt=paste(paste(set[-length(set)], res, collapse=" "), set[length(set)])
	return(txt)
}

bf=function(x, scale=c("lnL","logL")){
	list.obj=x 
#	list.obj must contain list of lnL (natural log) or logL (base 10) data
#   list.obj=list(a=c(-1,-4,-5), b=c(-10,-4,-5))
#	returns Bayes factor by computing the harmonic mean of likelihoods
	if("matrix"%in%class(list.obj)) list.obj=data.frame(list.obj)
	
	if(length(scale)>1) stop("'scale' of the data must be specified")
	
	if(any(sapply(list.obj, length)>1)) warning("Using harmonic mean estimator for marginal likelihood")
	hm=sapply(list.obj, hme, scale=scale)
	combn=combn(names(list.obj),2)
	t=data.frame(t(combn))
	if(scale=="lnL") {
		out=apply(combn, 2, function(x) {a=hm[[x[1]]]; b=hm[[x[2]]]; h=2*(hm[[x[1]]]-hm[[x[2]]]); return(c(a,b,h))})
		res=cbind(t, twice_lnBF=t(out))
		names(res)=c("m1", "m2", "lnL_m1", "lnL_m2", "twice_lnBF")
		interp=data.frame(twice_lnBF=c(2, 6, 10, Inf), interp=c("bare mention", "positive", "strong", "very strong"))
		rownames(interp)=c(0, 2, 6, 10)
		rr=res$twice_lnBF<0
		res$twice_lnBF=abs(res$twice_lnBF)

	} else if(scale=="logL") {
		out=apply(combn, 2, function(x) {a=hm[[x[1]]]; b=hm[[x[2]]]; h=hm[[x[1]]]-hm[[x[2]]]; return(c(a,b,h))})
		res=cbind(t, logBF=t(out))
		names(res)=c("m1", "m2", "logL_m1", "logL_m2", "logBF")
		interp=data.frame(logBF=c(1/2, 1, 2, Inf), interp=c("bare mention", "substantial", "strong", "decisive"))
		rownames(interp)=sprintf("%.1f",c(0, 1/2, 1, 2))
		rr=res$logBF<0
		res$logBF=abs(res$logBF)
	}
	rownames(res)=NULL
	resnum=res[,3:5]
	resfac=as.matrix(res[,1:2])
	ordres=matrix(NA, nrow=nrow(res), ncol=ncol(res))
	if(any(rr)){
		for(i in 1:nrow(res)){
			if(rr[i]) {
				resnum[i,]=resnum[i,c(2,1,3)]
				resfac[i,]=resfac[i,c(2,1)]
			} 
		}
	}
	res=cbind(resfac, resnum)
	res=res[order(res[,3], res[,4], decreasing=TRUE),]
	rownames(res)=NULL
	res$interp=sapply(1:nrow(res), function(i) as.character(interp$interp[min(which(abs(res[i,5])<interp[,1]))]))
	attr(interp, "reference")="Kass and Raftery 1995 J Amer Stat Assoc"
	attr(res, "interpretation")=interp
	attr(res, "txt")=.bf.interpreter(res)
	
	class(res)=c("bayesfactor", class(res))
	return(res)
}

print.bayesfactor=function(x, ...){
	cat("Bayes factor summary:\n\t")
	cat(attributes(x)$txt)
	cat("\n\n\n")
	cat(paste("Interpretation uses", sQuote(attributes(attributes(x)$interpretation)$reference)), ":\n", sep="")
	print(attributes(x)$interpretation)
	cat("\n\n")
	class(x)=class(x)[-which(class(x)=="bayesfactor")]
	print(x)
}

.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol




