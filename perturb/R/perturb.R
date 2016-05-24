perturb <- function(mod,pvars=NULL,prange=NULL,ptrans=NULL,pfac=NULL,uniform=FALSE,niter=100) {
	cutsp<-function(indx,tbl) {
		findInterval(runif(1),tbl[indx,],rightmost.closed=TRUE)
	}

	if (is.null(mod$call$formula)) stop("First argument does not contain a formula")
	stopifnot(is.list(pfac)||is.null(pfac))

	nms<-all.vars(terms(mod))
	stopifnot(all(pvars %in% nms))

	result<-NULL
	ncases<-length(get(nms[1]))
	frm<-deparse(formula(mod),width.cutoff = 500)
	result$formula<-frm
	allb<-coefficients(mod)

	# modify the formula
	if (length(pvars) > 0) {
		stopifnot(is.vector(get(pvars)))
		stopifnot(length(pvars)==length(prange))
		result$pvars<-pvars
		result$prange<-prange
		if (length(ptrans)>0) result$ptrans<-ptrans
		b<-make.names(c(nms,pvars),unique=TRUE)
		pvars.1<-b[(length(nms)+1):length(b)]
		for (i in 1:length(pvars)) {
			inp<-paste("\\<",pvars[i],"(\\>[^.]|$)",sep="")
			outp<-paste(pvars.1[i],"\\1",sep="")
			frm<-gsub(inp,outp,frm)
			ptrans<-gsub(inp,outp,ptrans)
		}
		result$ptrans2<-ptrans
	}
	if (length(pfac[[1]]) > 0) {
		rcls.tbl<-NULL
		pfac.1<-NULL
		if (is.list(pfac[[1]])) n<-length(pfac)
		else n<-1
		for (i in 1:n) {
			if (n == 1) lstnm<-pfac
			else lstnm<-pfac[[i]]
			stopifnot(all(lstnm[[1]] %in% nms))
			b<-make.names(c(nms,lstnm[[1]]),unique=TRUE)
			pfc<-b[(length(nms)+1):length(b)]
			inp<-paste("\\<",lstnm[[1]],"(\\>[^.]|$)",sep="")
			outp<-paste(pfc,"\\1",sep="")
			frm<-gsub(inp,outp,frm)
			rcls<-do.call("reclassify",lstnm)
			rcls.tbl<-c(rcls.tbl,list(rcls))
			pfac.1<-c(pfac.1,pfc)
		}
		result$reclassify.tables<-rcls.tbl
	}
	result$formula2<-frm
	mod$call$formula<-as.formula(frm)

	if (uniform) {
		ranexp<-"runif(ncases,-prange[i],prange[i])"
		result$distribution<-"uniform"
	}
	else {
		ranexp<-"rnorm(ncases,0,prange[i])"
		result$distribution<-"normal"
	}

	for (k in 1:niter) {
		# add random perturbances to pvars using values in prange
		if (length(prange)>0) {
			for (i in 1:length(prange)) {
				assign(pvars.1[i],get(pvars[i])+eval(parse(text=ranexp)))
			}
		}
		# re-execute the transformations
		for (trans in ptrans) eval(trans)
		# reclassify factors here
		if (length(pfac[[1]])>0) {
			for (i in 1:length(rcls.tbl)) {
				tbl<-rcls.tbl[[i]]$cum.reclass.prob
				tbl<-cbind(0,tbl)
				assign("tmpvar",as.numeric(get(rcls.tbl[[i]]$variable)))
				assign(pfac.1[i],sapply(tmpvar,cutsp,tbl))
				assign(pfac.1[i],as.factor(get(pfac.1[i])))
			}
		}
		# re-estimate the model using the perturbed variables
		mod2<-eval(mod$call)
		# collect the coefficients
		allb<-rbind(allb,coefficients(mod2))
	}
	# "allb" is the rowname value for the first row of allb; remove
	rownames(allb)<-NULL
	result$coeff.table<-allb
	class(result)<-"perturb"
	result
}

summary.perturb <-function(object,dec.places=3,full=FALSE,...) {
	coeffs<-object$coeff.table
	mysumm<-cbind(apply(coeffs,2,mean),apply(coeffs,2,sd),apply(coeffs,2,min),apply(coeffs,2,max))
	colnames(mysumm)<-c("mean","s.d.","min","max")
	object$coeff.table<-NULL
	object$summ<-mysumm
	object$dec.places<-dec.places
	object$full<-full
	dots<-substitute(expression(...))
	dots<-sub("^expression\\((.*)\\)$","\\1", deparse(dots))
	object$dots<-dots
	class(object)<-"summary.perturb"
	object
}

print.summary.perturb <-function(x,...) {
	if (x$full) {
		cat("formula:\n",x$formula,"\nformula2:\n",x$formula2,"\n\n")
	}
	if (length(x$pvars)>0) {
		cat("Perturb variables:\n")
		if (x$distribution=="uniform") {
			for (i in 1:length(x$pvars)) {
				prnt<-paste("uniform[",-round(x$prange[i],1),",",round(x$prange[i],1),"]",sep="")
				cat(x$pvars[i],"\t\t",prnt,"\n")
			}
		}
		else {
			for (i in 1:length(x$pvars)) {
				prnt<-paste("normal(0,",round(x$prange[i],1),")",sep="")
				cat(x$pvars[i],"\t\t",prnt,"\n")
			}
		}
		cat("\n")
	}
	if (length(x$ptrans)>0) {
		cat("Transformations:\n")
		for (trans in x$ptrans) {
			cat(trans,"\n")
		}
		if (x$full) {
			cat("\nTransformations2:\n")
			for (trans in x$ptrans2) {
				cat(trans,"\n")
			}
		}
		cat("\n")
	}
	if (!is.null(x$reclassify.tables)) {
		for (i in 1:length(x$reclassify.tables)) {
			if (x$dots=="") print(x$reclassify.tables[[i]],dec.places=x$dec.places,full=x$full,...)
			else eval(parse(text=paste("print(x$reclassify.tables[[i]],dec.places=x$dec.places,full=x$full,",x$dots,",...)")))
		}
	}
	cat("Impact of perturbations on coefficients:\n")
	#print(round(x$summ,x$dec.places),...)
	eval(parse(text=paste("print(round(x$summ,x$dec.places),",x$dots,",...)")))
}
