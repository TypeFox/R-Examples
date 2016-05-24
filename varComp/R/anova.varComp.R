anova.varComp <-
function (object, ..., test = c("KR", 'Satterthwaite'), L)
{
	cl=match.call(expand.dots=FALSE)
		diffLmat=function(qr0, X1)
		{
			z=(diag(1, nrow(X1))- tcrossprod(qr.Q(qr0)))%*%X1
			qrz=qr(zapsmall(z))    ### WARNING: this will depends on getOption("digits")
			 pivot <- qrz$pivot
			 oo <- order(pivot)
			 qr.R(qrz)[seq_len(qrz$rank), oo, drop=FALSE]
		}
	if('...'%in%names(cl)){
		if(!missing(L)) message("'L' will be ignored when '...' is given.")
		ddd.names=sapply(cl[['...']], deparse)
		ddd = list(...)
		nObj=length(ddd)
		nObj=nObj + 1L
		ddd[[nObj]] = object
		names(ddd) = c(ddd.names, deparse(substitute(object)))
		if(!all(sapply(ddd, inherits, what='varComp'))) stop("All `...` objects need to inherit the `varComp` class.")
		if(diff(range(sapply(ddd, nobs))) != 0 ) stop("The number of observations needs to be the same in all model fits.")
		n=nobs(ddd[[1L]])
		
		Xmats = lapply(ddd,  model.matrix, what='fixed')
		qrs = lapply(Xmats, qr)
		rks = sapply(qrs, '[[', 'rank')
		ork=order(rks)
		rks=rks[ork]
		
		ddd=ddd[ork]
		Xmats=Xmats[ork]
		qrs=qrs[ork]
		# (`raw.F`=drop(Fstat), `scale.F` = scale.overall, df1=rk, df2=F.ddf, `Pr(>F)`=drop(F.p))
		
		nas=rep(NA_real_, nObj)
		ans = data.frame(`F value`=nas, `Scale` = nas, numDF=nas, denDF=nas, `Pr(>F)` = nas, check.names=FALSE)
		rownames(ans) = names(ddd)
		
		if('model'%in%names(ddd[[1L]])){
			hasInt = attr(attr(ddd[[1L]]$model, 'terms'), 'intercept') == 1
		}else{
			ones = rep(1, n)
			hasInt = ( mean((tcrossprod(qr.Q(qrs[[1]]))%*%ones - ones)^2) < sqrt(.Machine$double.eps) )
		}
		if(hasInt && rks[1L] > 1L){
			Lmat = diffLmat(qr(rep(1,n)), Xmats[[1L]])
		}else{
			Lmat = diag(1, ncol(Xmats[[1L]]))
		}
		for(m in seq_len(nObj)){
			tmp = fixef(ddd[[m]], Lmat = Lmat, test=test)
			ans[m, ] = attr(attr(tmp, 'anova'), 'Overall')
			if(m < nObj)	Lmat = diffLmat(qrs[[m]], Xmats[[m+1L]])
		}
		return(ans)
	}
	
	### one object is given
	if(!missing(L)){
		tmp = fixef(object, Lmat=L, test=test)
	}else tmp = fixef(object, test=test)
	tmp.aov=attr(tmp, 'anova')
	ans = data.frame(`F value` = (tmp.aov[,'t value']*tmp.aov[,'Scale'])^2,
					 `Scale` = tmp.aov[,'Scale']^2, 
					 `numDF` = 1, 
					 `denDF` = tmp.aov[,'Df'], 
					 `Pr(>F)` = tmp.aov[, 'Pr(>|t|)'],
					 check.names=FALSE
					)
	if(missing(L)){
		message("Current implementation will test each fixed effect parameter separately when only one `varComp` object is provided.")
		X=model.matrix(object)
		n = nobs(object)
		if('model'%in%names(object)){
			hasInt = attr(attr(object$model, 'terms'), 'intercept') == 1
		}else{
			ones = rep(1, n)		
			hasInt = ( mean((tcrossprod(qr.Q(qr(X)))%*%ones - ones)^2) < sqrt(.Machine$double.eps) )
		}	
		if(hasInt){
			tmp = fixef(object, Lmat = diffLmat(qr(rep(1,n)), X), test=test)
			tmp.aov=attr(tmp, 'anova')
		}
	}
	ans = rbind(ans, attr(tmp.aov, 'Overall'))
	if(!missing(L) && !is.null(rownames(L))) rownames(ans) = c(rownames(L), 'Overall')
	return(ans)

}
