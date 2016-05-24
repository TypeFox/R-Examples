reclassify <- function(varname,pcnt=NULL,adjust=TRUE,bestmod=TRUE,min.val=.1,diag=NULL,unif=NULL,dist=NULL,assoc=NULL) {
	# subroutines used below:
	exptab <- function(pcnt) {
		result<-NULL
		mdl2<-NULL
		if (length(pcnt)==1) {
			init.pcnt<-(100-pcnt)/(100*(ncat-1))
			init.pcnt<-matrix(init.pcnt,ncat,ncat)
			diag(init.pcnt)<-pcnt/100
			mdl<-"diag1"
		}
		else if (length(pcnt)==ncat) {
			init.pcnt<-(100-pcnt)/(100*(ncat-1))
			init.pcnt<-replicate(ncat,cbind(init.pcnt))
			diag(init.pcnt)<-pcnt/100
			mdl<-"diag"
		}
		else if (length(pcnt)==ncat^2) {
			init.pcnt<-as.matrix(pcnt)
			dim(init.pcnt)<-c(ncat,ncat)
			init.pcnt<-prop.table(init.pcnt,1)
			mdl<-"diag1+unif"
			mdl2<-"diag1+dist"
		}
		else stop("Invalid argument for pcnt option")
		init.tbl<-diag(freq)%*%init.pcnt
		init.tbl[init.tbl < min.val]<-min.val
		result$init.pcnt<-init.pcnt
		result$init.tbl<-init.tbl
		result$mdl<-mdl
		result$mdl2<-mdl2
		result
	}

	best.model <- function(init.tbl,ncat,mdl,mdl2) {
		result<-NULL
		ncase<-length(init.tbl)
		pm<-init.tbl
		dim(pm)<-c(ncase,1)
		orig<-gl(ncat,ncat,ncase)
		dest<-gl(ncat,1,ncase)
		diag<-as.factor((as.numeric(orig)==as.numeric(dest))*(as.numeric(orig)))
		diag1<-as.numeric(as.numeric(orig)==as.numeric(dest))
		dist<-abs(as.numeric(orig)-as.numeric(dest))
		unif<-as.numeric(orig)*as.numeric(dest)
		m<-suppressWarnings(glm(as.formula(paste("pm~orig+dest+",mdl)),family=poisson(),maxit=1000))
		if (!is.null(mdl2)) {
			m2<-suppressWarnings(glm(as.formula(paste("pm~orig+dest+",mdl2)),family=poisson(),maxit=1000))
			if (m2$deviance < m$deviance ) {
				m<-m2
				mdl<-mdl2
			}
		}
		cf.all<-coef(m)
		mmat<-model.matrix(as.formula(paste("~",mdl)))
		cf.names<-colnames(mmat)[-1]
		cf<-cf.all[cf.names]
		mmat<-as.matrix(mmat[,-1])
		paras<-mmat%*%cf
		result$model<-mdl
		result$coefs<-cf
		result$paras<-paras
		result
	}

	assocPat <- function(freq,ncat,diag=NULL,unif=NULL,dist=NULL,assoc=NULL) {
		result<-NULL
		cf <- NULL
		tbl<-diag(freq)

		ncase<-length(tbl)
		dim(tbl)<-c(ncase,1)
		paras<-rep(0,ncase)
		orig<-gl(ncat,ncat,ncase)
		dest<-gl(ncat,1,ncase)
		rmat<-model.matrix(~orig)
		rmat<-rmat[,attr(rmat,"assign")==1]
		cmat<-model.matrix(~dest)
		cmat<-cmat[,attr(cmat,"assign")==1]
		eqmain<-rmat+cmat

		if (!is.null(assoc)) {
			dim(assoc)<-NULL
			paras<-assoc
		}
		else {
			if (!is.null(diag)) {
				paras <- paras + as.numeric(orig==dest)*diag
				cf <- cbind(diag,cf)
				colnames(cf)[1]<-"diag"
			}
			if (!is.null(unif)) {
				paras <- paras + as.numeric(orig)*as.numeric(dest)*unif
				cf <- cbind(unif,cf)
				colnames(cf)[1]<-"unif"
			}
			if (!is.null(dist)) {
				paras <- paras + abs(as.numeric(orig)-as.numeric(dest))*dist
				cf <- cbind(dist,cf)
				colnames(cf)[1]<-"dist"
			}
		}
		ll<-suppressWarnings(glm(tbl~eqmain+offset(paras),family=poisson(),maxit=1000))
		pred<-ll$fitted.values
		dim(pred)<-c(ncat,ncat)
		result$pred<-pred
		result$iter<-ll$iter
		result$converged<-ll$converged
		result$coefs<-cf
		result
	}

	### reclassify starts here
	if (all(is.null(pcnt),is.null(assoc),is.null(diag),is.null(unif),is.null(dist))) {
		stop("Either pcnt=, assoc= or one of diag=, unif=, or dist= must be specified")
	}
	if (!adjust) bestmod<-FALSE
	result<-NULL
	if (is.character(varname)) {
		result$variable<-varname
		varname<-get(varname)
	}
	else {
		result$variable<-substitute(varname)
	}
	stopifnot(is.factor(varname))
	freq<-table(varname)
	ncat<-length(freq)

	if (!is.null(pcnt)) {
		if (is.character(pcnt)) {pcnt<-get(pcnt)}
		etab<-exptab(pcnt)
		result$exptab<-etab
		if (bestmod) {
			bm<-best.model(etab$init.tbl,ncat,etab$mdl,etab$mdl2)
			result$bestmod<-bm
			pred<-assocPat(freq,ncat,assoc=bm$paras)
		}
		else if (adjust) {
			init.tbl<-etab$init.tbl
			init.tbl[init.tbl==0]<-1e-12
			init.tbl<-log(init.tbl)
			init.tbl<-(init.tbl+t(init.tbl))/2
			pred<-assocPat(freq,ncat,assoc=init.tbl)
		}
		else {
			rcd<-prop.table(etab$init.tbl,1)
			result$reclass.prob<-rcd
			rcdcum<-t(apply(rcd,1,cumsum))
			colnames(rcdcum)<-levels(varname)
			result$cum.reclass.prob<-rcdcum
			class(result)<-"reclassify"
			return(result)
		}
	}
	else if (!is.null(assoc)){
		if (is.character(assoc)) {assoc<-get(assoc)}
		stopifnot(length(assoc)==ncat^2)
		dim(assoc)<-c(ncat,ncat)
		assoc<-(assoc+t(assoc))/2
		rownames(assoc)<-levels(varname)
		colnames(assoc)<-levels(varname)
		pred<-assocPat(freq,ncat,assoc=assoc)
		result$assoc<-assoc
	}
	else if (!is.null(diag) | !is.null(unif) | !is.null(dist)) {
		pred<-assocPat(freq,ncat,diag=diag,unif=unif,dist=dist)
		result$coefs<-pred$coefs
	}

	if (!pred$converged) {
		stop(cat("The requested association pattern could not be fitted to",
		substitute(varname),"in",pred$iter,"iterations."))
	}
	rownames(pred$pred)<-levels(varname)
	colnames(pred$pred)<-levels(varname)
	rcd<-prop.table(pred$pred,1)
	result$fitted.table<-pred$pred

	result$reclass.prob<-rcd
	rcdcum<-t(apply(rcd,1,cumsum))
	colnames(rcdcum)<-levels(varname)
	result$cum.reclass.prob<-rcdcum
	class(result)<-"reclassify"
	result
}

print.reclassify <-function(x,dec.places=3,full=FALSE,...) {
	prnt.tbl<-function(obj,header,lbl=TRUE,...){
		cat(header,"\n")
		if (lbl) cat("original\treclassified\n")
		print(round(obj,digits=dec.places),quote=FALSE,...)
		cat("\n")
	}

	## print.reclassify starts here
	cat("\nVariable",x$variable,"to be reclassified ")
	if (full) {
		if (!is.null(x$exptab)) {
			prnt.tbl(x$exptab$init.pcnt,"using the following initial probabilities:")
			prnt.tbl(addmargins(x$exptab$init.tbl),"Initial expected table based on these probabilities:")
		}
		if (!is.null(x$bestmod)) {
			prnt.tbl(x$bestmod$coefs,"Best model for the intial expected table:",lbl=FALSE)
		}
		else if (!is.null(x$assoc)) {
			prnt.tbl(x$assoc,"using the following log pattern of association:")
		}
		else if (!is.null(x$coefs)) {
			prnt.tbl(x$coefs,"using the following model coefficients:",lbl=FALSE)
		}
		if (!is.null(x$fitted.table)) {
			prnt.tbl(addmargins(x$fitted.table),"Table generated with specified pattern of association:")
		}
	}
	cat("Reclassification probabilities:\n")
	cat("original\treclassified\n")
	print(round(x$reclass.prob,digits=dec.places),quote=FALSE,...)
	cat("\n")
	if (full) {
		cat("Cumulative reclassification probabilities:\n")
		cat("original\treclassified\n")
		print(round(x$cum.reclass.prob,digits=dec.places),quote=FALSE,...)
		cat("\n")
	}
}
