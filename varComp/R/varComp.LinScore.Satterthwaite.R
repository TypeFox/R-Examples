varComp.LinScore.Satterthwaite <-
function( # null.fit,
	all.scores, infoMat, null, w, tr1, n, ...)
{
#Actual Satterthwaite approximation for LinScore. Argument definitions are the same as in varComp.LinScore.SSAS155, except the LIkLI, tau.idx, and non.pd are not used. 
  nK=length(all.scores)
  rqforms=(all.scores+tr1/2)  # n/2 * ratio of qforms
#  
  nonNull=seq_len(nK)[-null]
#  
  Phi=infoMat[null, nonNull, drop=FALSE]  
  Delta=infoMat[null, null, drop=FALSE]   
#  
  var.score=infoMat[nonNull, nonNull, drop=FALSE] - crossprod(Phi, solve(Delta, Phi))
  mean.score=tr1[nonNull]*.5 - crossprod(Phi, solve(Delta, tr1[null]*.5-rqforms[null]))
#  
  sum.score=sum(w*mean.score)
  var.sum=drop(crossprod(w, var.score%*%w))
#  
  scale=var.sum/2/sum.score
  df=2*sum.score^2/var.sum
#
  obs.stat=sum(rqforms[-null]*w)
  pval=pchisq(obs.stat/scale, df, lower.tail=FALSE)
#  
  ans=list(statistic=c(`sum of ratio of quadratic forms`=obs.stat), 
           p.value=pval, 
           alternative='greater', 
           parameter=c(scale=scale, df=df), 
           null.value=structure(numeric(length(nonNull)), names=sprintf('variance component %d', nonNull)), #null.fit=null.fit, 
           method='Satterthwaite Approximation of Profiled Variance Component Test'
          )
  class(ans)='htest'
  ans  
}
