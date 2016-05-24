varComp.LinScore.test <-
function(control, n, tau.idx, LIkLI, tr1, infoMat, all.scores, non.pd)  #, ...)
{
#### external ref: null.fit, X, Y, K, non.pd

#Compute weights for LinScore; do test when the null has no other variance components; pass to varComp.LinScore.approximation when the null has other variance components. 
#i.	control: the control object from varCompTest.control.
#ii.	null: the same as in varComp.test. 
#iii.	n: residual sample size.
#v.	tau.idx: the complement of null. 
#vi.	LIkLI: the same as in varComp.VM03.test.
#vii.	tr1: the same as varComp.VM03.test
#viii.	obs.score: observed score vector for the alternative\null
#ix.	infoMat: the same as varComp.VM03.test
#x.	all.scores: all scores
#xi.	...: place holder for unused arguments. 
  null = seq_along(all.scores)[-tau.idx]
  nNull=length(null)
  wts=control$wt
  obs.score = all.scores[tau.idx]
  if(length(null)==0L) {
    methods=control$method[[1]] 
    if(!any(methods%in%c('AS155', 'exact', 'Davies', 'SSAS155')))
      message('LinScore test method was changed to the exact method "AS155" when null=integer(0L). ')
    methods='AS155'
  }else{
    methods=control$method[[2]]
    tmp=intersect(methods, c('AS155', 'exact', 'Davies'))
    if(length(tmp)>0L){
      warning(sprintf('LinScore test method(s) [%s] was changed to "SSAS155" when null!=integer(0L).', paste(tmp,collapse=', ')))
      methods=unique(c('SSAS155', setdiff(methods, tmp)))
    }    
  }
  n.tau=length(tau.idx); nK=length(LIkLI)
  ans=vector('list', length(wts)); names(ans)=wts
  
  if(nNull>0L){
    partNegHess=infoMat[tau.idx,tau.idx,drop=FALSE]-
      infoMat[tau.idx,-tau.idx,drop=FALSE]%*%solve(infoMat[-tau.idx,-tau.idx,drop=FALSE],infoMat[-tau.idx,tau.idx,drop=FALSE])
    invNegHess=tryCatch(solve(partNegHess), error=function(e){partNegHess<<-as.matrix(nearPD(partNegHess)$mat);solve(partNegHess)})
  }else{
    partNegHess=infoMat; invNegHess=solve(infoMat)
  }
  if(any(grepl('^Wd', wts))){ ## Wald variant
    inv.info=solve(infoMat)
    mle=inv.info%*%all.scores
  }
    
  for(w.i in seq_along(wts)){
    tmpwts=gsub('^Wd','', wts[w.i])
    if(tmpwts!=wts[w.i]){  ## Wald
      bak.partNegHess=partNegHess
      partNegHess=inv.info[tau.idx, tau.idx, drop=FALSE]
    }
    w=switch(tmpwts,
        EqWt=rep(1/n.tau, n.tau),
        InvSTD=local({tmp=1/sqrt(diag(as.matrix(partNegHess))); tmp/sum(tmp)}),
        InvSqrtV=local({tmp=pmax(0, base::colSums(t(backsolve(chol(partNegHess), diag(1,nK))))); tmp/sum(tmp)}),
        MinVar=solve.QP(Dmat=partNegHess, dvec=numeric(n.tau), Amat=cbind(1, diag(1,n.tau)), bvec=rep(1:0,c(1,n.tau)), meq=1, factorized=FALSE)$solution 
      )
    if(tmpwts!=wts[w.i]){  ## Wald
      w=drop(w%*%inv.info[tau.idx, , drop=FALSE])
      partNegHess=bak.partNegHess
      lin.form=sum(w*all.scores)
    }else
      lin.form=sum(w*obs.score)

    this.ans=vector('list', length(methods)); names(this.ans)=methods
    for(m.i in seq_along(methods)){
      if(nNull==0L){
        stopifnot(methods[m.i]=='AS155')
        #pval=davies(0, eigen(n*Reduce('+', mapply('*',w,LIkLI[tau.idx],SIMPLIFY=FALSE)),TRUE,TRUE)$val - sum(tr1[tau.idx]*w)-lin.form*2, acc=control$acc, lim=control$lim) ## when length(null)==0L, tau.idx subsetting is not needed. 
        pval=davies(0, eigen(n*Reduce('+', mapply('*',w,LIkLI,SIMPLIFY=FALSE)),TRUE,TRUE)$val - sum(tr1*w)-lin.form*2, acc=control$acc, lim=control$lim)
        this.ans[[m.i]]=structure(list(statistic=c(`linear form`=lin.form), 
                                       p.value=pval$Qq, 
                                       parameter=c(w=w,
                                                   acc=control$acc, lim=control$lim,
                                                   trace=pval$trace, 
                                                   ifault=pval$ifault
                                                      ), 
                                       method='Exact variance component test', 
                                       alternative='greater',
                                       null.value=structure(numeric(n.tau), names=sprintf('variance component %d', tau.idx))
                                      ),
                                  class='htest'
                                 )
      }else{
        this.ans[[m.i]]=do.call(
			what = paste0('varComp.LinScore.', methods[m.i]),
			args=list(null=null, 
				w=w, lin.form=lin.form, LIkLI=LIkLI, tr1=tr1
				,infoMat=infoMat, all.scores=all.scores, n=n, tau.idx=tau.idx, non.pd=non.pd, control=control
				)
			)
      }
      ans[[w.i]]=this.ans
    }
  }
  class(ans)='varComp.LinScore.test'
  ans
}
