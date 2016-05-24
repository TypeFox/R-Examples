"MLE.beta"<-function(X.list, ped, beta=NULL, nUSdam=NULL, nUSsire=NULL){

  off_records<-match(as.character(X.list$id[as.numeric(names(X.list$X))]), as.character(ped[,1]))
  dam<-match(ped[,2][off_records], X.list$id)
  sire<-match(ped[,3][off_records], X.list$id)

  if(sum(is.na(dam)==FALSE)==0 & sum(is.na(sire)==FALSE)==0){
     warning("no pedigree structure - cannot estimate beta")
     C <- diag(rep(1,length(unique(X.list$beta_map))))
     B.start <- rep(0,length(unique(X.list$beta_map)))
  }else{

    nUS<-c(0,0)

    if(length(nUSdam)>0){nUS[1]<-nUSdam[1]}
    if(length(nUSsire)>0){nUS[2]<-nUSsire[1]}

    nbeta<-c(ncol(X.list$X[[1]]$XDus), ncol(X.list$X[[1]]$XDs), ncol(X.list$X[[1]]$XSus), ncol(X.list$X[[1]]$XSs), ncol(X.list$X[[1]]$XDSus), ncol(X.list$X[[1]]$XDSs))

    if(length(beta)==0){
      beta<-matrix(0, length(unique(X.list$beta_map)),1)
    }else{
      beta<-as.matrix(beta)
    }

    nind<-length(X.list$id)
    ndam<-unlist(lapply(X.list$X, function(x){length(x$dam.id)}))
    nsire<-unlist(lapply(X.list$X, function(x){length(x$sire.id)}))
    dam_pos<-1:length(X.list$X)
    sire_pos<-1:length(X.list$X)
    par_pos<-1:length(X.list$X)

    X<-lapply(X.list$X, function(x){list(DS=NULL, S=NULL, D=NULL)})
    mergeN<-lapply(X.list$X, function(x){x$mergeN})

    for(i in 1:length(X.list$X)){

      if(is.null(X.list$merge)==FALSE){
        for(m in 1:length(X.list$merge)){
          if(X.list$merge[m]<=sum(nbeta[1:2])){
            if(X.list$mergeUS[m]==1){
              mergeN[[i]][,m][1]<-mergeN[[i]][,m][1]+nUS[1]
            }
            if(X.list$mergeUS[m]==2){
              mergeN[[i]][,m][2]<-mergeN[[i]][,m][2]+nUS[1]
            }
          }else{ 
            if(X.list$mergeUS[m]==1){
              mergeN[[i]][,m][1]<-mergeN[[i]][,m][1]+nUS[2]
            }
            if(X.list$mergeUS[m]==2){
              mergeN[[i]][,m][2]<-mergeN[[i]][,m][2]+nUS[2]
            }
          }
        }
      }

      if(is.na(dam[i])){
        if(X.list$X[[i]]$restdam.id[length(X.list$X[[i]]$restdam.id)]>nind){
          dam_pos[i]<-length(X.list$X[[i]]$restdam.id)
        }else{
          dam_pos[i]<-NA
        }
      }else{
        dam_pos[i]<-match(dam[i], X.list$X[[i]]$dam.id)
      }

      if(is.na(sire[i])){
       if(X.list$X[[i]]$restsire.id[length(X.list$X[[i]]$restsire.id)]>nind){
          sire_pos[i]<-length(X.list$X[[i]]$restsire.id)
        }else{
          sire_pos[i]<-NA
        }
      }else{
        sire_pos[i]<-match(sire[i], X.list$X[[i]]$sire.id)
      }

      par_pos[i]<-(nsire[i]*(dam_pos[i]-1))+sire_pos[i]

      if(sum(nbeta[1:2])>0){
        X[[i]]$D<-matrix(NA, ndam[i], sum(nbeta[1:2]))
      }
      if(sum(nbeta[3:4])>0){
        X[[i]]$S<-matrix(NA, nsire[i], sum(nbeta[3:4]))
      }
      if(sum(nbeta[5:6])>0){
        X[[i]]$DS<-matrix(NA, ndam[i]*nsire[i], sum(nbeta[1:6])) 
      }
      if(nbeta[1]>0){
        X[[i]]$D[,1:nbeta[1]]<-X.list$X[[i]]$XDus
      }
      if(nbeta[2]>0){
        X[[i]]$D[,nbeta[1]+(1:nbeta[2])]<-X.list$X[[i]]$XDs
      }
      if(nbeta[3]>0){
        X[[i]]$S[,1:nbeta[3]]<-X.list$X[[i]]$XSus
      }
      if(nbeta[4]>0){
        X[[i]]$S[,nbeta[3]+(1:nbeta[4])]<-X.list$X[[i]]$XSs
      }
      if(nbeta[5]>0){
        X[[i]]$DS[,sum(nbeta[1:4])+(1:nbeta[5])]<-X.list$X[[i]]$XDSus
      }
      if(nbeta[6]>0){
        X[[i]]$DS[,sum(nbeta[1:5])+(1:nbeta[6])]<-X.list$X[[i]]$XDSs
      }
      if(sum(nbeta[5:6])>0){
        if(sum(nbeta[1:2])>0){
          X[[i]]$DS[,1:sum(nbeta[1:2])]<-apply(X[[i]]$D, 2, rep, each=nsire[i])
          X[[i]]$D<-NULL
        }
        if(sum(nbeta[3:4])>0){
          X[[i]]$DS[,sum(nbeta[1:2])+(1:sum(nbeta[3:4]))]<-apply(X[[i]]$S, 2, rep, ndam[i])
          X[[i]]$S<-NULL
        }
      }
    }

    optim.out <- optim(beta, beta.loglik, method = "BFGS",control = list(fnscale = -1), hessian = TRUE, X=X, dam_pos=dam_pos, sire_pos=sire_pos, par_pos=par_pos, beta_map=X.list$beta_map, merge=X.list$merge, mergeN=mergeN, nUS=nUS)
 
    C <- solve(-1 * optim.out$hessian)
    B.start <- optim.out$par
  }
list(beta=B.start, C=C)
} 


