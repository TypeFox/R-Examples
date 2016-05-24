fixef.varComp <-
function(object, Lmat, alpha=.05, test=c('KR', 'Satterthwaite'), ...)
{
## S3 method for reporting fixed effect parameters from class varComp. 
# object: varComp object
# Lmat: The same as in satterth, default to all parameters (i.e., identity Lmat)
# alpha: Type I error level requested. 
	if(is.null(test[1L])) return(coef(object, 'fixed'))
	test=match.arg(test)

	if('Y'%in%names(object)){
		Y = object$Y
	}else if ('model'%in%names(object)){
		Y = model.response(object$model)
	}else stop("response variable is not recored.")

  
  X=model.matrix(object, what='fixed')
  K=model.matrix(object, what='K')[object$parms>0]
  if(missing(Lmat)) {Lmat=diag(1, ncol(X)); rownames(Lmat)=colnames(X)}
  if(!is.matrix(Lmat))	Lmat = matrix(Lmat, ncol=length(Lmat), dimnames = list("", names(Lmat)))
  
  if(!is.null(colnames(Lmat)) && !is.null(colnames(X))){
	if(ncol(Lmat) < ncol(X)){
		L=matrix(0, nrow(Lmat), ncol(X))
		colnames(L) = colnames(X)
		rownames(L) = rownames(Lmat)
		L[, colnames(Lmat)]=Lmat
		Lmat=L
	}else if (ncol(Lmat) == ncol(X)){
		Lmat = Lmat[, colnames(X), drop=FALSE]
	}else stop("`Lmat` has more columns than the number of fixed effect parameters.")	
  }
  if(ncol(Lmat) != ncol(X)) stop("`Lmat` has incorrect number of columns.")
  if(is.null(rownames(Lmat))) rownames(Lmat) = rep('', nrow(Lmat))
  
  this.V=vcov(object, what='Y')
  this.Vbet=if(ncol(X)>0L) solve(crossprod(X, solve(this.V,X))) else matrix(NA_real_, 0L, 0L)
  this.bet=this.Vbet%*%crossprod(X, solve(this.V, Y))
  
  est=drop(Lmat%*%this.bet)
  LVL=tcrossprod(Lmat%*%this.Vbet, Lmat)
  svd.LVL=if(length(LVL)>0L) svd(LVL) else list(d=numeric(0L), u=matrix(NA_real_, 0L,0L), v=matrix(NA_real_,0L,0L))
  rk=sum(svd.LVL$d>sqrt(.Machine$double.eps))
    if(test=='Satterthwaite'){
	  
	  ses=sqrt(diag(LVL))
	  t.dfs=numeric(nrow(Lmat))
	  for(i in seq_len(nrow(Lmat))) 
		t.dfs[i]=satterth.varComp(object=object, Lmat[i,,drop=FALSE], Vbet=this.Vbet, 
						  #svd.VLbet=list(d=svd.LVL$d[i], u=svd.LVL$u[i,,drop=FALSE], v=svd.LVL$v[i,,drop=FALSE]),
						  X=X, K=K, V=this.V )
	#  t.dfs=apply(Lmat, 1L, satterth, object=object, Vbet=this.Vbet, X, K=K, V=this.V)  ## X has naming conflict
		scaleF=rep(1, nrow(Lmat)); scale.overall = if(nrow(Lmat)>0L) 1 else numeric(0L)
	  #debug(satterth.varComp)
	  LVLI=tcrossprod(sweep(svd.LVL$v, 2, ifelse(svd.LVL$d>sqrt(.Machine$double.eps), 1/svd.LVL$d, 0), '*'), svd.LVL$u) 
	  Fstat=if(nrow(Lmat)>0L) crossprod(est, LVLI%*%est)/rk  else numeric(0L)
	  F.ddf=satterth.varComp(object, Lmat=Lmat, Vbet=this.Vbet, svd.VLbet=svd.LVL, K=K, V=this.V, X)
	  F.p=pf(Fstat, rk, F.ddf, lower.tail=FALSE)
	}else if (test=='KR') {
		scaleF = ses = t.dfs =numeric(nrow(Lmat))
		for(i in seq_len(nrow(Lmat))){
			tmp = KR.varComp(object=object, Lmat=Lmat[i,,drop=FALSE], Vbet=this.Vbet, X = X, K = K, V = this.V)
			t.dfs[i]=tmp[[1L]]
			ses[i] = sqrt(diag(Lmat[i,,drop=FALSE]%*%attr(tmp, 'vcov.beta')%*%t(Lmat[i,,drop=FALSE])))
			scaleF[i] = attr(tmp, 'Scale')
		}
		tmp=KR.varComp(object=object, Lmat=Lmat, Vbet=this.Vbet, svdVLbet = svd.LVL, X = X, K = K, V = this.V)
		scale.overall = attr(tmp, 'Scale')
		Fstat = attr(tmp, 'F value')
		F.ddf = if(length(tmp)>0L) tmp[[1L]] else numeric(0L)
		F.ndf = attr(tmp, 'numDF')
		if(length(F.ndf)>0L && F.ndf != rk) warning("numerical instability detected on the rank of `Lmat`")
		F.p = attr(tmp, 'Pr(>F)')
	}
	ll = est - ses/sqrt(scaleF) * qt(1-alpha/2, t.dfs) 
	ul = est + ses/sqrt(scaleF) * qt(1-alpha/2, t.dfs) 
	tstats = est / ses 
	t.pval = pt(abs(tstats)*sqrt(scaleF), t.dfs, lower.tail=FALSE)*2
	

  ans = est
  attrs = cbind(`Std. Error`=ses, Lower=ll, Upper=ul, `t value`=tstats*sqrt(scaleF), `Scale` = sqrt(scaleF), Df=t.dfs, `Pr(>|t|)`=t.pval)
  attr(attrs, 'Overall')=cbind(`F value`=drop(Fstat*scale.overall), `Scale` = scale.overall, numDF=if(nrow(Lmat)>0L)rk else integer(0L), denDF=F.ddf, `Pr(>F)`=drop(F.p))
  names(est) = rownames(attrs)=rownames(Lmat)
  rownames(attr(attrs, 'Overall')) = if(nrow(Lmat)>0L) 'Overall' else character(0L)
  attr(ans, 'anova') = attrs
  class(ans) = 'varCompFixEf'
  ans
}

print.varCompFixEf = function(x,...)
{
	cat("Individual fixef effect estimates:\n")
	tmp=cbind(Estimate=x, attr(x, 'anova'))
	print(tmp)
	cat("\nOverall fixed effect contrast:\n")
	print(attr(attr(x, 'anova'), 'Overall'))
	cat('\n')
	invisible(x)
}