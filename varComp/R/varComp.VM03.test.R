varComp.VM03.test <-
function(control, infoMat, tau.idx, LIkLI, tr1, n, LIy, all.scores)# , ...)
{#	varComp.VM03.test: Actual testing function for VM03 test
#i.	control: the control object from varCompTest.control. 
#ii.	infoMat: Negative expected Hessian matrix.
#iii.	tau.idx: Integer vector that indexes the additional variance components of the alternative hypothesis compared to the null. That is, this is the complement of the null argument of varComp.test.
#iv.	LIkLI: A list of matrices being  sqrt(V^{-1})' K_i sqrt(V ^{-1})
#v.	tr1: Vector of trace terms. 
#vi.	n: residual sample size. 
#vii.	LIy: sqrt(V^{-1})'y
#viii.	null: the same as in varComp.test. 
#ix.	all.scores: vector of all score statistics. 
#x.	...: place holder for unused arguments. 
  null = seq_along(all.scores)[-tau.idx]
  nNull=length(null)
  n.tau=length(tau.idx)
  if(length(null)==0L) methods=control$method[[1]] else methods=control$method[[2]]
  if(nNull>0L && length(setdiff(methods, 'pboot'))==0L) { ## FIXME: implement pboot
 	warning("parametric bootstrap has not been implemented yet for the situation when the null model is not a linear model")
   null.htest=structure(list(p.value=NA_real_),class='htest')
    ans=vector('list', length(methods)); names(ans)=methods
    for(m in methods) ans[[m]]=null.htest
    class(ans)='varComp.VM03.test'; return(ans) 
  }
  #sumLIkLI=mapply('*',w,LIkLI[tau.idx])
#  invNegHess=solve(infoMat)[tau.idx, tau.idx,drop=FALSE]  ## this should be the same as below
  if(nNull>0L){ #warning('VM03 is not seriously tested when null!=integer(0) and should be used with extreme caution!')
    Phi=infoMat[-tau.idx, tau.idx, drop=FALSE]
    Delta=infoMat[-tau.idx, -tau.idx, drop=FALSE]
    partNegHess=infoMat[tau.idx,tau.idx,drop=FALSE]- crossprod(Phi, solve(Delta, Phi))
    invNegHess=tryCatch(solve(partNegHess), error=function(e){partNegHess<<-as.matrix(nearPD(partNegHess)$mat);solve(partNegHess)})
    
    mean.nonNull.score=drop(crossprod(Phi, solve(Delta, all.scores[-tau.idx])))

    adj.nonNull.score=all.scores[tau.idx] - mean.nonNull.score
    Amat=diag(1,n.tau); bvec=rep(0,n.tau)
    fact=1; i=0L
    repeat{
      i=i+1
      tmp=try( solve.QP(Dmat=invNegHess/fact, dvec=drop(invNegHess%*%all.scores[tau.idx])/fact, Amat=Amat, bvec=bvec, meq=0L, factorized=FALSE)$value )
      if(class(tmp)!='try-error'){
        qp.term=tmp*fact
        break
      }else if(i==20L){
        null.htest=structure(list(p.value=NA_real_),class='htest')
        ans=vector('list', length(methods)); names(ans)=methods
        for(m in methods) ans[[m]]=null.htest
        class(ans)='varComp.VM03.test'; 
        warning('Quadratic solver error')
        return(ans) 
      }else fact=fact*2
    }
    SSstat0=max(0, -2*qp.term      )
  }else{
    partNegHess=infoMat; invNegHess=solve(infoMat)
    Amat=diag(1,n.tau); bvec=rep(0,n.tau)
    SSstat=function(z){
      S=.5* (sapply(LIkLI, function(lik) drop(crossprod(z, lik%*%z)/crossprod(z)))*n-tr1 )
      SInvNegHess=crossprod(invNegHess,S)
      #max(0, -2*solve.QP(Dmat=invNegHess, dvec=drop(SInvNegHess), Amat=Amat, bvec=bvec, meq=0L, factorized=FALSE)$value )
      sol=drop(solve.QP(Dmat=invNegHess, dvec=drop(SInvNegHess), Amat=Amat, bvec=bvec, meq=0L, factorized=FALSE)$solution)
      max(0, crossprod(sol, invNegHess%*%sol))
    }
    SSstat0=drop(SSstat(LIy))
    nsim=control$nsim
  }
  ans=vector('list', length(methods)); names(ans)=methods
  for(m.i in seq_along(methods)){
    if(methods[m.i]=='pboot'){
      if(nNull==0L){
        seed=get.seed()
        SSstats=replicate(nsim-1L, SSstat(rnorm(n)))
        pval=(1+sum(SSstats>=SSstat0))/(nsim)
        ans[[m.i]]=structure(list(statistic=c(VM03=SSstat0), p.value=pval, null.stats=SSstats, seed=seed, parameters=c(nsim=nsim)), class='htest')      
      }else{
        .NotYetImplemented()
      }
    }else if(methods[m.i]=='ChiBarSq'){
      if(nNull==0L){
        seed=get.seed()
        pval=pchibarsq(SSstat0, partNegHess, lower.tail=FALSE)
        ans[[m.i]]=structure(list(statistic=c(VM03=SSstat0), p.value=pval, seed=seed, parameters=partNegHess), class='htest')
      }else{
        if(length(tau.idx)==1L){  ## this is a misnomer, as the mixture is not a chi-bar-square
          pval=if(SSstat0==0) 1 else pnorm(sqrt(SSstat0), mean.nonNull.score*sqrt(invNegHess), lower.tail=FALSE)
          ans[[m.i]]=structure(list(statistic=c(VM03=SSstat0), p.value=pval, parameters=c(mean=mean.nonNull.score, var=partNegHess)), class='htest')
        }else{
          .NotYetImplemented()
        }
      }
    }else if(methods[m.i]=='Bound'){
      if(nNull==0L){ ## the same as ChiBarSq
        seed=get.seed()
        pval=pchibarsq(SSstat0, partNegHess, lower.tail=FALSE)
        ans[[m.i]]=structure(list(statistic=SSstat0, p.value=pval, seed=seed, parameters=partNegHess), class='htest')
      }else{
        .NotYetImplemented()
      }
    }else if(methods[m.i]=='SSChiBarSq'){
      if(nNull==0L){
        seed=get.seed()
        SSstats=replicate(nsim, SSstat(rnorm(n)))
        m=mean(SSstats); s=sd(SSstats)
        theo.mom=mchibarsq(partNegHess, order=1:2)
        theo.sd=sqrt(theo.mom[2]-theo.mom[1]^2)
        SSstat0.adj=(SSstat0-m)/s*theo.sd+theo.mom[1]
        pval=pchibarsq(SSstat0.adj, partNegHess, lower.tail=FALSE)
        ans[[m.i]]=structure(list(statistic=c(VM03=SSstat0, VM03.adj=SSstat0.adj), p.value=pval, 
                                  null.stats=SSstats, seed=seed, 
                                  parameters=c(nsim=nsim, emp.mean=m, emp.sd=s, 
                                                theo.mean=theo.mom[1], theo.sd=theo.sd)), 
                                  
                             class='htest')    
      }else{
        .NotYetImplemented()
        seed=get.seed()
        pval=pchibarsq(SSstat0, partNegHess, lower.tail=FALSE)
        ans[[m.i]]=structure(list(statistic=SSstat0, p.value=pval, seed=seed, parameters=partNegHess), class='htest')
      }      
    }else stop('Undefined methods for VM03')
  }
  class(ans)='varComp.VM03.test'
  ans
}
