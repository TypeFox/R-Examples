varComp.SS95.test <-
function(control,infoMat, tau.idx, LIkLI, LIy, tr1, n, all.scores)# , ...)
{
#Actual testing function for SS95 test. Arguments are the same as in varComp.VM03.test. 
	null = seq_along(all.scores)[-tau.idx]
  nNull=length(null); n.tau=length(tau.idx)
  if(length(null)==0L) methods=control$method[[1]] else methods=control$method[[2]]
  if(nNull>0L && length(setdiff(methods, 'pboot'))==0L) { ## FIXME: implement pboot
	warning("parametric bootstrap has not been implemented yet for the situation when the null model is not a linear model")
    null.htest=structure(list(p.value=NA_real_),class='htest')
    ans=vector('list', length(methods)); names(ans)=methods
    for(m in methods) ans[[m]]=null.htest
    class(ans)='varComp.SS95.test'; return(ans) 
  }
  #sumLIkLI=mapply('*',w,LIkLI[tau.idx])
  if(nNull>0L){
    Phi=infoMat[-tau.idx, tau.idx, drop=FALSE]
    Delta=infoMat[-tau.idx, -tau.idx, drop=FALSE]
    partNegHess=infoMat[tau.idx,tau.idx,drop=FALSE]-
      infoMat[tau.idx,-tau.idx,drop=FALSE]%*%solve(infoMat[-tau.idx,-tau.idx,drop=FALSE],infoMat[-tau.idx,tau.idx,drop=FALSE])
    invNegHess=tryCatch(solve(partNegHess), error=function(e){partNegHess<<-as.matrix(nearPD(partNegHess)$mat);solve(partNegHess)})
    mean.nonNull.score=drop(crossprod(Phi, solve(Delta, all.scores[-tau.idx])))
    adj.nonNull.score=all.scores[tau.idx] - mean.nonNull.score
    SS95stat0=max(0, -2*solve.QP(partNegHess, drop(adj.nonNull.score), diag(1,n.tau), rep(0,n.tau))$value)
  }else{
    partNegHess=infoMat; invNegHess=solve(infoMat)
    Amat=diag(1, n.tau); bvec=numeric(n.tau)
    SS95stat=function(z){
      S=.5* (sapply(LIkLI[tau.idx], function(lik) drop(crossprod(z, lik%*%z)/crossprod(z)))*n-tr1[tau.idx] )
      tvec=solve.QP(partNegHess, drop(S), Amat, bvec)$solution
      max(0, crossprod(tvec, partNegHess%*%tvec))
    }
    SS95stat0=drop(SS95stat(LIy))
  }
  nsim=control$nsim
  ans=vector('list', length(methods)); names(ans)=methods
  for(m.i in seq_along(methods)){
    if(methods[m.i]=='pboot'){
      if(nNull==0L){
        seed=get.seed()
        SS95stats=replicate(nsim-1L, SS95stat(rnorm(n)))
        pval=(1+sum(SS95stats>=SS95stat0))/(nsim)
        ans[[m.i]]=structure(list(statistic=c(SS95=SS95stat0), p.value=pval, null.stats=SS95stats, seed=seed, parameters=c( nsim=nsim)), class='htest')
      }else{
        .NotYetImplemented()
      }
    }else if(methods[m.i]=='ChiBarSq'){
      seed=get.seed()
      pval=pchibarsq(SS95stat0, invNegHess, lower.tail=FALSE)
      ans[[m.i]]=structure(list(statistic=c(SS95=SS95stat0), p.value=pval, seed=seed, parameters=invNegHess), class='htest')
    }else stop('Undefined methods for SS95')
  }
  class(ans)='varComp.SS95.test'
  ans
}
