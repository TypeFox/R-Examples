sfem <- function(Y,K=2:6,model='AkjB',method='reg',crit='icl',maxit=50,eps=1e-6,init='kmeans',nstart=25,Tinit=c(),kernel='',disp=F,l1=0.1,l2=0,nbit=2){
  call = match.call()
  if (max(l1)>1 || min(l1)<0 || max(l2)>1 || min(l2)<0) stop("Parameters l1 and l2 must be within [0,1]\n",call.=FALSE)
  else{
    bic = aic = icl = c()
    RES = list()
    try(res0 <- fem(Y,K,init=init,nstart=nstart,maxit=maxit,eps=eps,Tinit=Tinit,model=model,kernel=kernel,method='reg',crit=crit))
    for (i in 1:length(l1)){
      try(RES[[i]] <- fem.sparse(Y,res0$K,model=res0$model,maxit=15,eps=eps,Tinit=res0$P,l1=l1[i],l2=l2,nbit=nbit))
      try(bic[i] <- RES[[i]]$bic)
      try(aic[i] <- RES[[i]]$aic)
      try(icl[i] <- RES[[i]]$icl)
      #try(fish[i] <- RES[[i]]$fish)
    }
    if (crit=='bic'){ id_max = which.max(bic); crit_max = RES[[id_max]]$bic}
    if (crit=='aic'){ id_max = which.max(aic); crit_max = RES[[id_max]]$aic}
    if (crit=='icl'){ id_max = which.max(icl); crit_max = RES[[id_max]]$icl}
    #if (crit=='fisher'){ id_max = which.max(diff(fish)); crit_max = RES[[id_max]]$fish}
    res = RES[[id_max]]
    res$call = res0$call
    res$plot = res0$plot
    res$plot$l1 = list(aic=aic,bic=bic,icl=icl,l1=l1)
    res$call = call
    res$crit = crit
    res$l1 = l1[id_max]
    res$l2 = l2
    if (disp){ if (length(l1)>1) cat('The best sparse model is with lambda =',l1[id_max],'(',crit,'=',crit_max,')\n')
    else cat('The sparse model has a lambda =',l1[id_max],'(',crit,'=',crit_max,')\n')}
  }
  class(res)='fem'
  res
}
