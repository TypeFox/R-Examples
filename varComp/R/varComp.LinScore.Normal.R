varComp.LinScore.Normal <-
function(	# null.fit, 
	all.scores, lin.form, infoMat, null, w, tr1, n, LIkLI, tau.idx,  ...)
{
#Actual normal approximation for LinScore. Argument definitions are the same as in varComp.LinScore.SSAS155, except the non.pd is not used. 
  nK=length(all.scores)
  nonNull=seq_len(nK)[-null]
  Phi=infoMat[null, nonNull, drop=FALSE]  
  Delta=infoMat[null, null, drop=FALSE]   

  var.nonNull.score=infoMat[nonNull, nonNull, drop=FALSE] - crossprod(Phi, solve(Delta, Phi))
  mean.nonNull.score= crossprod(Phi, solve(Delta, all.scores[null]))
  mean.obs.score=sum(mean.nonNull.score*w)
  var.obs.score=drop(crossprod(w, var.nonNull.score%*%w))
  
  mean.theo.score=0 
  var.theo.score=drop(crossprod(w, infoMat[nonNull, nonNull, drop=FALSE]%*%w))
  
  adjusted.obs.score=sqrt(1/var.obs.score)*(lin.form-mean.obs.score)
  
  pval=pnorm(adjusted.obs.score, lower.tail=FALSE)
  pval0=pnorm(lin.form/sqrt(var.obs.score), lower.tail=FALSE)
  structure(list(statistic0=c(z=lin.form/sqrt(var.obs.score)), 
                 p.value0=pval0, 
                 statistic=c(z=adjusted.obs.score), 
                 p.value=pval,
                 parameter=c(lin.form=lin.form, sd=sqrt(var.obs.score)), 
                 method='Profiled Variance Component Test Based on Asymptotic Normality', 
                 alternative='greater',
                 null.value=structure(numeric(length(nonNull)), names=sprintf('variance component %d', tau.idx)) 	#, null.fit=null.fit
                ),
            class='htest'
           )
}
