"reordXlist"<-function(X.list, marker.type="MSW"){

    nind<-length(X.list$id)
      
    if(TRUE%in%(X.list$X[[1]]$dam.id>nind)){
      USdam<-TRUE
    }else{
      USdam<-FALSE
    }   

    if(TRUE%in%(X.list$X[[1]]$sire.id>nind)){
      USsire<-TRUE
    }else{
      USsire<-FALSE
    }

  npar<-c(ncol(X.list$X[[1]]$XDus), ncol(X.list$X[[1]]$XDs), ncol(X.list$X[[1]]$XSus), ncol(X.list$X[[1]]$XSs), ncol(X.list$X[[1]]$XDSus), ncol(X.list$X[[1]]$XDSs))

  for(off in 1:length(X.list$X)){

    if(is.null(X.list$X[[off]]$mmS)){
      Sorder<-rep(1, length(X.list$X[[off]]$sire.id))
      Dorder<-rep(1, length(X.list$X[[off]]$dam.id))
    }else{
      if(marker.type!="AFLP"){
        Sorder<-X.list$X[[off]]$mmS
        Dorder<-X.list$X[[off]]$mmD
      }else{
        if(is.null(X.list$X[[off]]$G)){
          Sorder<-rep(1, length(X.list$X[[off]]$sire.id))
          Dorder<-rep(1, length(X.list$X[[off]]$dam.id))
        }else{
          Sorder<-matrix(-X.list$X[[off]]$G, length(X.list$X[[off]]$sire.id), length(X.list$X[[off]]$dam.id))
          Dorder<-colSums(Sorder)
          Sorder<-rowSums(Sorder)
        }
      }
    }

    if(USdam==TRUE){
       Dorder[which(X.list$X[[off]]$dam.id>nind)]<-999
    }

    if(USsire==TRUE){
       Sorder[which(X.list$X[[off]]$sire.id>nind)]<-999
    }
    
    if(length(X.list$X[[off]]$restdam.id)!=length(X.list$X[[off]]$dam.id)){
       Dorder[which(X.list$X[[off]]$dam.id%in%X.list$X[[off]]$restdam.id==FALSE)]<-9999
    }

     if(length(X.list$X[[off]]$restsire.id)!=length(X.list$X[[off]]$sire.id)){
       Sorder[which(X.list$X[[off]]$sire.id%in%X.list$X[[off]]$restsire.id==FALSE)]<-9999
    }

    restDorder<-order(Dorder[match(X.list$X[[off]]$restdam.id, X.list$X[[off]]$dam.id)])
    restSorder<-order(Sorder[match(X.list$X[[off]]$restsire.id, X.list$X[[off]]$sire.id)])


    Dorder<-order(Dorder)
    Sorder<-order(Sorder)

    X.list$X[[off]]$mmD<-X.list$X[[off]]$mmD[Dorder]
    X.list$X[[off]]$mmS<-X.list$X[[off]]$mmS[Sorder]

    X.list$X[[off]]$dam.id<-X.list$X[[off]]$dam.id[Dorder]
    X.list$X[[off]]$sire.id<-X.list$X[[off]]$sire.id[Sorder]

    X.list$X[[off]]$restdam.id<-X.list$X[[off]]$restdam.id[restDorder]
    X.list$X[[off]]$restsire.id<-X.list$X[[off]]$restsire.id[restSorder]

    if(npar[1]>0){
      X.list$X[[off]]$XDus<-matrix(X.list$X[[off]]$XDus[Dorder,], ncol=npar[1], dimnames=dimnames(X.list$X[[off]]$XDus))
    }
    if(npar[2]>0){
      X.list$X[[off]]$XDs<-matrix(X.list$X[[off]]$XDs[Dorder,],  ncol=npar[2], dimnames=dimnames(X.list$X[[off]]$XDs))
    }
    if(npar[3]>0){
      X.list$X[[off]]$XSus<-matrix(X.list$X[[off]]$XSus[Sorder,],  ncol=npar[3], dimnames=dimnames(X.list$X[[off]]$XSus))
    }
    if(npar[4]>0){
      X.list$X[[off]]$XSs<-matrix(X.list$X[[off]]$XSs[Sorder,],  ncol=npar[4], dimnames=dimnames(X.list$X[[off]]$XSs))
    }
    if(npar[5]>0){
      X.list$X[[off]]$XDSus<-matrix(X.list$X[[off]]$XDSus[c(t(outer(length(X.list$X[[off]]$sire.id)*(Dorder-1), Sorder, "+"))),],  ncol=npar[5], dimnames=dimnames(X.list$X[[off]]$XDSus))
    }
    if(npar[6]>0){
      X.list$X[[off]]$XDSs<-matrix(X.list$X[[off]]$XDSs[c(t(outer(length(X.list$X[[off]]$sire.id)*(Dorder-1), Sorder, "+"))),],  ncol=npar[6], dimnames=dimnames(X.list$X[[off]]$XDSs))
    }
    if(is.null(X.list$X[[off]]$G)==FALSE){
      X.list$X[[off]]$G<-as.matrix(X.list$X[[off]]$G[c(t(outer(length(X.list$X[[off]]$sire.id)*(Dorder-1), Sorder, "+")))])
    }
  }  
X.list
}

