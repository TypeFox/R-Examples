varComp.RLRT.test <-
function(control,negHess, alt.fit, null.fit, tau.idx)  #, k,null,...)
{
#Actual testing function for RLRT. 
#i.	control. the control object from varCompTest.control. 
#ii.	negHess: Negative expected Hessian matrix. 
#iii.	alt.fit: Object from varComp under alternative hypothesis. 
#iv.	null.fit: Object form varComp under null hypothesis. 
#v.	tau.idx: complement of null. 
#vi.	k: list of corrected correlation matrices. 
#vii.	null: the same as the null in varComp.test. 
#viii.	...: place holder for unused arguments. 
	
  nNull=length(null.fit$random.labels); 
  null = seq_along(alt.fit$random.labels)[-tau.idx]
  n.tau=length(tau.idx)
  k = alt.fit$working.cor
  if(n.tau==1L) methods=control$method[[1]] else methods=control$method[[2]]
#  if(nNull>0L) { ## FIXME: implement pboot
#    null.htest=structure(list(p.value=NA_real_),class='htest')
#    ans=vector('list', length(methods)); names(ans)=methods
#    for(m in methods) ans[[m]]=null.htest
#    class(ans)='varComp.RLRT.test'; return(ans) 
#  }
  #require(RLRsim)
  LRstat0=max(0, 2*(alt.fit$PREML-null.fit$PREML))
  
  nsim=control$nsim
  ans=vector('list', length(methods)); names(ans)=methods
  for(m.i in seq_along(methods)){
    if(methods[m.i]=='pboot' && nNull==0L && n.tau==1L)  methods[m.i]='exact' ## C & R (2004) situation
    if(methods[m.i]=='pboot'){
      .NotYetImplemented()
    }else if( methods[m.i]=='exact'){
      if(n.tau!=1L){
        #warning('RLRT is not implemented for multiple variance components')
        ans[[m.i]]=structure(list(statistic=c(RLR=LRstat0), p.value=NA_real_, parameters=c(nsim=nsim), null.stats=NA_real_), class='htest')
      }else{
        K.half=suppressWarnings(t(chol(k[[tau.idx]], pivot=TRUE)))
        if(nNull>0L)  warning('RLRT when the null model has other variance component than the error variance is only approximate!')
        
        oo=order(attr(K.half, 'pivot'))                                       ## fixed by LQ on 050412
        K.des=K.half[oo, seq_len(attr(K.half,'rank'))]
        LRstats=RLRTSim(matrix(NA_real_,nrow(K.des),0L), Z=K.des, sqrt.Sigma=diag(1,ncol(K.des)), nsim=nsim-1L)
        pval=(1+sum(LRstats>=LRstat0))/(nsim)
        ans[[m.i]]=structure(list(statistic=c(RLR=LRstat0), p.value=pval, parameters=c(nsim=nsim),null.stats=LRstats
                                  ,null.fit=null.fit, alt.fit=alt.fit
                                  ), class='htest')
      }
    }else stop('Unrecognized RLRT method')
  }
  class(ans)='varComp.RLRT.test'
  ans
}
