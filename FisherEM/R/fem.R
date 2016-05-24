fem <- function(Y,K=2:6,model='AkjB',method='reg',crit='icl',maxit=50,eps=1e-6,init='kmeans',nstart=25,Tinit=c(),kernel='',disp=F){
  call = match.call()
  MOD = c('DkBk','DkB','DBk','DB','AkjBk','AkjB','AkBk','AkB','AjBk','AjB','ABk','AB','all')
  MET = c('svd','reg','gs')
  KER = c('','sigmoid','linear','rbf')
  CRI = c('bic','aic','icl')
  INIT = c('user','random','kmeans','hclust')
  if (any(!model%in%MOD)) stop("Invalid model name\n",call.=FALSE)
  if (any(!method%in%MET)) stop("Invalid method name\n",call.=FALSE)  
  if (any(!kernel%in%KER)) stop("Invalid kernel name\n",call.=FALSE)
  if (any(!crit%in%CRI)) stop("Invalid criterion name\n",call.=FALSE)
  if (any(!init%in%INIT)) stop("Invalid initialization name\n",call.=FALSE)
  if (init=='hclust' & nrow(Y)>5000) stop('Too much data for this initialization',,call.=FALSE)
  # if (K>=ncol(Y)) stop("K must be strictly less than the number of variables",call.=FALSE)
  if (nrow(Y)<=ncol(Y) & method=='gs') stop("n<<p case: use method REG or SVD ",call.=FALSE)
  if (length(model)==1) if (model=='all') model = c('DkBk','DkB','DBk','DB','AkjBk','AkjB','AkBk','AkB','AjBk','AjB','ABk','AB')
  resultat = vector("list",length=length(K)*length(model))
    bic = aic = icl = rep(NA,length(K)*length(model))
    for (k in 1:length(K)){
      if (disp) cat('>> K =',K[k],'\n')
      for(i in 1:length(model)){
        try(resultat[[(k-1)*length(model)+i]]<-fem.main(Y,K[k],init=init,nstart=nstart,maxit=maxit,eps=eps,Tinit=Tinit,model=model[i],kernel=kernel,method=method))
        if (length(resultat[[(k-1)*length(model)+i]])>0){
          try(bic[(k-1)*length(model)+i]<-resultat[[(k-1)*length(model)+i]]$bic) 
          try(aic[(k-1)*length(model)+i]<-resultat[[(k-1)*length(model)+i]]$aic)
          try(icl[(k-1)*length(model)+i]<-resultat[[(k-1)*length(model)+i]]$icl)
	  #try(fish[(k-1)*length(model)+i]<-resultat[[(k-1)*length(model)+i]]$fish)
	  if (disp){
          	if (crit=='bic') cat(MOD[i],'\t:\t bic =',resultat[[(k-1)*length(model)+i]]$bic,'\n')
          	if (crit=='aic') cat(MOD[i],'\t:\t aic =',resultat[[(k-1)*length(model)+i]]$aic,'\n')
          	if (crit=='icl') cat(MOD[i],'\t:\t icl =',resultat[[(k-1)*length(model)+i]]$icl,'\n')
	  	#if (crit=='fisher') cat(MOD[i],'\t:\t fish =',resultat[[(k-1)*length(model)+i]]$fish)
	}
      }}
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
    res$plot = list(bic=matrix(bic,ncol=nm,byrow=T),aic=matrix(aic,ncol=nm,byrow=T),icl=matrix(icl,ncol=nm,byrow=T),K=K,models=model)
    res$call = call
    class(res)='fem'
    res
}

