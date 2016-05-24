varComp.LinScore.SSAS155 <-
function(	#null.fit, 
	all.scores, lin.form, infoMat, null, w, tr1, n, LIkLI, tau.idx, non.pd, control, ...)
{
#Actual shifted and scaled AS155
#i.	null.fit: object from varComp under the null. 
#ii.	all.scores: All scores.
#iii.	lin.form: test statistic. 
#iv.	infoMat: information matrix. 
#v.	null: the same as in varComp.test.
#vi.	w: vector of weights. 
#vii.	tr1: vector of traces.
#viii.	n: residual sample size
#ix.	LIkLI: list of matrices, each being sqrt(V^{-1})' K_i sqrt(V ^{-1})
#x.	tau.idx: complement of null.
#xi.	non.pd: A logical indicating that the infoMat is not positive definite. This only applies when the information argument of varComp.test is not `EI`. 
#control: Optional control object, primarily for getting extremely small p-values
#xii.	...: place holder for unused arguments.
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
  
  adjusted.obs.score=sqrt(var.theo.score/var.obs.score)*(lin.form-mean.obs.score) + mean.theo.score
  
  pval=davies(0, eigen(n*Reduce('+', mapply('*',w,LIkLI[tau.idx],SIMPLIFY=FALSE)),TRUE,TRUE)$val - sum(tr1[tau.idx]*w)-adjusted.obs.score*2, acc=control$acc, lim=control$lim)
  pval0=davies(0, eigen(n*Reduce('+', mapply('*',w,LIkLI[tau.idx],SIMPLIFY=FALSE)),TRUE,TRUE)$val - sum(tr1[tau.idx]*w)-lin.form*2,acc=control$acc, lim=control$lim)
  structure(list(statistic0=c(`linear form`=lin.form), 
                 p.value0=pval0$Qq, 
                 statistic=c(`adjusted linear form`=adjusted.obs.score), 
                 p.value=pval$Qq,
                 parameter=c(w=w, 
                                 acc=control$acc, 
                                 lim=control$lim, 
                                 ifault=pval$ifault, 
                                 trace=pval$trace, 
                                 non.pd=non.pd
                                ), 
                 method='Profiled Variance Component Test through Shifted & Scaled AS155', 
                 alternative='greater',
                 null.value=structure(numeric(length(nonNull)), names=sprintf('variance component %d', tau.idx)) #, null.fit=null.fit
                ),
            class='htest'
           )
}
