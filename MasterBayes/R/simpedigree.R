simpedigree<-function(PdP, beta=NULL, nUS=NULL){
 
  USdam<-PdP$USdam
  USsire<-PdP$USsire

  PdP$USdam<-FALSE
  PdP$USsire<-FALSE

  X.list<-getXlist(PdP)

  unique_id<-as.character(X.list$id)

  ped<-cbind(as.character(unique_id), rep(NA, length(unique_id)), rep(NA, length(unique_id)))
  
  nd<-length(X.list$X[[1]]$XDs)/length(X.list$X[[1]]$dam.id)
  ns<-length(X.list$X[[1]]$XSs)/length(X.list$X[[1]]$sire.id)
  nds<-length(X.list$X[[1]]$XDSs)/(length(X.list$X[[1]]$dam.id)*length(X.list$X[[1]]$sire.id))
  
  nrecords<-length(PdP$id)
           
  nbeta<-length(unique(X.list$beta_map))                   # number of parameters 
  if(X.list$beta_map[1]==-999){nbeta<-0}

  if(length(beta)==0){
    beta<-matrix(0, nbeta,1)
  }else{
    if(length(beta)==nbeta){
      beta<-as.matrix(beta)
    }else{
      warning("length beta does not equal the number of parameters in the model")
    }
  }

  if(length(USdam)==1){
    if(USdam[1]==FALSE){
      USdam<-NULL
    }else{
      USdam<-rep(1,nrecords)
    }
  }

  if(length(USsire)==1){
    if(USsire[1]==FALSE){
      USsire<-NULL
    }else{
      USsire<-rep(1,nrecords)
    }
  }

  nusd<-length(unique(USdam))            # number of unsampled female catgeories
  nuss<-length(unique(USsire))           # number of unsampled male catgeories

  X<-lapply(X.list$X, function(x){matrix(NA, length(x$dam.id)*length(x$sire.id), length(X.list$beta_map))})

  for(i in 1:length(X.list$X)){
    if(nd>0){
        X[[i]][,1:nd]<-apply(X.list$X[[i]]$XDs, 2, rep, each=length(X.list$X[[i]]$sire.id))
    }
    if(ns>0){
        X[[i]][,nd+1:ns]<-apply(X.list$X[[i]]$XSs, 2, rep, length(X.list$X[[i]]$dam.id))
    }
    if(nds>0){
        X[[i]][,(nd+ns)+1:nds]<-X.list$X[[i]]$XDSs
    }
  }


    if(length(beta)==0){
      alpha<-lapply(X.list$X, function(x){rep(1, length(x$dam.id)*length(x$sire.id))})
    }else{
      if(is.null(X.list$merge)){
        alpha<-lapply(X, function(x){exp(x%*%beta[X.list$beta_map])})
      }else{
        alpha<-list()
        for(i in 1:length(X.list$X)){
           beta_tmp<-beta[X.list$beta_map]
           for(m in 1:length(X.list$merge)){
             beta_tmp[m]<-inv.logit(beta_tmp[m])
             n1<-X.list$X[[i]]$mergeN[,m][1]  
             n2<-X.list$X[[i]]$mergeN[,m][2] 
             beta_tmp[m]<-(beta_tmp[m]/n1)/((beta_tmp[m]/n1)+((1-beta_tmp[m])/n2))
             beta_tmp[m]<-logit(beta_tmp[m])
           }
           alpha[[i]]<-exp(X[[i]]%*%beta_tmp)
        }
      }       
    }
 
  dam<-1:length(X.list$X)
  sire<-1:length(X.list$X)

  patmat<-lapply(alpha, function(x){as.numeric(sample(1:length(x),1,prob=x))})


  for(i in 1:length(X.list$X)){ 
    dam[i]<-ceiling(patmat[[i]]/length(X.list$X[[i]]$sire.id))
    sire[i]<-patmat[[i]]-((as.numeric(dam[i])-1)*(length(X.list$X[[i]]$sire.id)))
    dam[i]<-unique_id[X.list$X[[i]]$dam.id[as.numeric(dam[i])]]
    sire[i]<-unique_id[X.list$X[[i]]$sire.id[as.numeric(sire[i])]]
  }

  ped[,2][as.numeric(names(X.list$X))]<-dam
  ped[,3][as.numeric(names(X.list$X))]<-sire

  output<-list(ped=ped, USsire.data=NULL, USsire.formula=NULL, USdam.data=NULL, USdam.formula=NULL)

  if(nusd>0){
    output$USdam.data<-as.matrix(rep(0,nrecords))
    colnames(output$USdam.data)[1]<-"USdam.data"
    output$USdam.formula<-expression(varPed("USdam.data", gender="Female", relational=FALSE, restrict=0))
    for(us_cat in 1:nusd){
      d_delete<-sample(subset(1:nrecords, USdam==unique(USdam)[us_cat] &  PdP$sex=="Female" &  PdP$offspring==0), USdam[us_cat])  
      offspring_in_cat<-match(PdP$id[which(PdP$offspring==1 & USdam==unique(USdam)[us_cat])], ped[,1])
      output$USdam.data[d_delete]<-1
      offspring_in_cat<-offspring_in_cat[which(ped[,2][offspring_in_cat]%in%PdP$id[d_delete]==TRUE)]
      output$ped[,2][offspring_in_cat]<-NA
    }
  }

  if(nuss>0){
    output$USsire.data<-as.matrix(rep(0,nrecords))
    colnames(output$USsire.data)[1]<-"USsire.data"
    output$USsire.formula<-expression(varPed("USsire.data", gender="Male", relational=FALSE, restrict=0))
    for(us_cat in 1:nuss){
      s_delete<-sample(subset(1:nrecords, USsire==unique(USsire)[us_cat] &  PdP$sex=="Male" &  PdP$offspring==0), USsire[us_cat])  
      offspring_in_cat<-match(PdP$id[which(PdP$offspring==1 & USsire==unique(USsire)[us_cat])], ped[,1])
      output$USsire.data[s_delete]<-1
      offspring_in_cat<-offspring_in_cat[which(ped[,3][offspring_in_cat]%in%PdP$id[s_delete]==TRUE)]
      output$ped[,3][offspring_in_cat]<-NA      
    }
  }


output
}
