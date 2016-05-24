funFEM <-
function(fd,K=2:6,model='AkjBk',crit='bic',init='hclust',Tinit=c(),maxit=50,eps=1e-6,disp=FALSE,lambda=0,graph=FALSE){
  call = match.call()
  MOD = c('DkBk','DkB','DBk','DB','AkjBk','AkjB','AkBk','AkB','AjBk','AjB','ABk','AB','all')
  CRI = c('bic','aic','icl')
  INIT = c('user','random','kmeans','hclust')
  if (any(!model%in%MOD)) stop("Invalid model name\n",call.=FALSE)
  if (any(!crit%in%CRI)) stop("Invalid criterion name\n",call.=FALSE)
  if (any(!init%in%INIT)) stop("Invalid initialization name\n",call.=FALSE)
  # if (init=='hclust' & nrow(Y)>5000) stop('Too much data for this initialization',,call.=FALSE)
  # if (K>=ncol(Y)) stop("K must be strictly less than the number of variables",call.=FALSE)
  if (length(model)==1) if (model=='all') model = c('DkBk','DkB','DBk','DB','AkjBk','AkjB','AkBk','AkB','AjBk','AjB','ABk','AB')
  resultat = vector("list",length=length(K)*length(model))
  bic = aic = icl = ll= nbprm = rep(NA,length(K)*length(model))
  it = 1
  for (k in 1:length(K)){
    if (disp) cat('>> K =',K[k],'\n')
    for(i in 1:length(model)){
      try(resultat[[it]]<-.FunFEM.main(fd,K[k],init=init,maxit=maxit,eps=eps,Tinit=Tinit,model=model[i],lambda=lambda,graph=graph))
      if (length(resultat[[it]])>0){
        try(bic[it]<-resultat[[it]]$bic) 
        try(aic[it]<-resultat[[it]]$aic)
        try(icl[it]<-resultat[[it]]$icl)
        try(nbprm[it]<-resultat[[it]]$nbprm)
        try(ll[it]<-resultat[[it]]$ll)
        #try(fish[(k-1)*length(model)+i]<-resultat[[(k-1)*length(model)+i]]$fish)
        if (disp){
          if (crit=='bic') cat(model[i],'\t:\t bic =',resultat[[it]]$bic,'\n')
          if (crit=='aic') cat(model[i],'\t:\t aic =',resultat[[it]]$aic,'\n')
          if (crit=='icl') cat(model[i],'\t:\t icl =',resultat[[it]]$icl,'\n')
          #if (crit=='fisher') cat(MOD[i],'\t:\t fish =',resultat[[(k-1)*length(model)+i]]$fish)
        }
        it = it + 1
      }
    }
    if (disp) cat('\n')
  }
  if (crit=='bic'){ id_max = which.max(bic); crit_max = resultat[[id_max]]$bic}
  if (crit=='aic'){ id_max = which.max(aic); crit_max = resultat[[id_max]]$aic}
  if (crit=='icl'){ id_max = which.max(icl); crit_max = resultat[[id_max]]$icl}
  #if (crit=='fisher'){ id_max = which.max(diff(fish)); crit_max = resultat[[id_max]]$fish}
  res = resultat[[id_max]]
  if (disp) cat('The best model is',res$model,'with K =',res$K,'(',crit,'=',crit_max,')\n')
  res$crit = crit
  nm = length(model)
  res$plot = list(bic=matrix(bic,ncol=nm,byrow=T),aic=matrix(aic,ncol=nm,byrow=T),icl=matrix(icl,ncol=nm,byrow=T),K=K,nbprm=matrix(nbprm,ncol=nm,byrow=T),ll=matrix(ll,ncol=nm,byrow=T),models=model)
  res$call = call
  class(res)='fem'
  res
}
