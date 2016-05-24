pfda <-
function(Y,cls,model='AkjBk',crit='bic',cv.fold=10,kernel='',display=FALSE){
  MOD = c('DkBk','DkB','DBk','DB','AkjBk','AkjB','AkBk','AkB','AjBk','AjB','ABk','AB','all')
  KER = c('','sigmoid','linear','rbf')
  CRI = c('bic','cv')
  if (any(!model%in%MOD)) stop("Invalid model name\n",call.=FALSE)
  if (any(!kernel%in%KER)) stop("Invalid kernel name\n",call.=FALSE)
  if (any(!crit%in%CRI)) stop("Invalid criterion name\n",call.=FALSE)
  
  if (length(model)==1) if (model=='all') model = c('DkBk','DkB','DBk','DB','AkjBk','AkjB','AkBk','AkB','AjBk','AjB','ABk','AB')
  resultat = vector("list",length=length(model))
  
  if (crit == "cv"){
    CRIT = matrix(NA,cv.fold,length(model))
    for (v in 1:cv.fold){
      test = sample(nrow(Y),round(nrow(Y)/4))
      for(i in 1:length(model)){
        prms = .pfda.main(Y[-test,],cls[-test],model=model[i])
        res = predict(prms,Y[test,])
        CRIT[v,i] = sum(res$cl==cls[test])/length(cls[test])
      }
    }
    CRIT = colMeans(CRIT)
    if (display) {names(CRIT) = model; print(CRIT)}
  }
  
  if (crit == "bic"){
    CRIT = rep(NA,length(model))
    for(i in 1:length(model)){
      CRIT[i] = .pfda.main(Y,cls,model=model[i])$bic
    }
    if (display) {names(CRIT) = model; print(CRIT)}
  }
  
  if (display) cat('* Selected model:',model[which.max(CRIT)],'\n')
  res = .pfda.main(Y,cls,model=model[which.max(CRIT)])
}
