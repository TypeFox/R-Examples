"MLE.popsize"<-function(X.list, USdam=FALSE, USsire=FALSE, nUS=NULL, ped=NULL){
 

  if(length(USdam)==1){
    if(USdam==TRUE){
      nbetaD<-1
      USdam<-rep(1, length(X.list$X))
    }else{
      nbetaD<-0
    }
  }else{  
    if(length(USdam)!=length(X.list$X)){stop("length of USdam does not equal number of offspring")}
    nbetaD<-length(unique(USdam))
  }

  if(USsire[1]=="USdam"){
    USsiredam<-TRUE
    USsire<-USdam
  }else{  
    USsiredam<-FALSE
  }

  if(length(USsire)==1){
    if(USsire==TRUE){
      nbetaS<-1
      USsire<-rep(1, length(X.list$X))
    }else{
      nbetaS<-0
    }
  }else{
    if(length(USsire)!=length(X.list$X)){stop("length of USsire does not equal number of offspring")}
    nbetaS<-length(unique(USsire))
  }

  if(length(nUS)==0){
    nUS<-matrix(1E-5, (nbetaD+nbetaS*(USsiredam==FALSE)),1)
  }else{
    if(length(nUS)!=(nbetaD+nbetaS*(USsiredam==FALSE))){
      warning("beta is wrong size in popsize.loglik")
      stop()
    }else{
      nUS<-as.matrix(nUS)
    }
  }

    lower<-rep(1E-5,nbetaD+nbetaS*(USsiredam==FALSE))
    upper<-rep(1E+7,nbetaD+nbetaS*(USsiredam==FALSE))

    nsire<-unlist(lapply(X.list$X, function(x){length(x$sire.id)}))
    ndam<-unlist(lapply(X.list$X, function(x){length(x$dam.id)}))
    nudam<-unlist(lapply(X.list$X, function(x){length(x$restdam.id)}))
    nusire<-unlist(lapply(X.list$X, function(x){length(x$restsire.id)}))

    npar<-ndam*nsire 

    if(is.null(ped)==FALSE){
      if(sum(is.na(ped[,2])==FALSE)==0 & nbetaD>0){
        ped<-NULL
      }
    }

    if(is.null(ped)==FALSE){
      if(sum(is.na(ped[,3])==FALSE)==0 & nbetaS>0){
        ped<-NULL
      }
    }

    if(is.null(ped)==FALSE){
      ped<-ped[match(X.list$id[as.numeric(names(X.list$X))], ped[,1]),]
    }

    X<-lapply(X.list$X, function(x){list(N=NULL,G=NULL)})

        for(i in 1:length(X)){
          nsampD<-ndam[i]-(nbetaD>0)
          nsampS<-nsire[i]-(nbetaS>0)

          sampD<-c(1:npar[i])
          sampS<-c(1:npar[i])

          if(nbetaD>0){
            sampD<-c(1:npar[i])[-c((1:nsire[i])+((nudam[i]-1)*nsire[i]))]
          }
          if(nbetaS>0){
            sampS<-c(1:npar[i])[-c(((0:(ndam[i]-1))*nsire[i])+nusire[i])]
          }
          DandS<-intersect(sampS, sampD)
          DnotS<-setdiff(sampD, sampS)
          SnotD<-setdiff(sampS, sampD)
          notDS<-(nbetaD>0)*(nbetaS>0)*(nusire[i]+((nudam[i]-1)*nsire[i]))
          X[[i]]$N<-c(nsampD*nsampS, nsampS,  nsampD, 1)   
          if(is.null(ped)){
            X[[i]]$G<-c(sum(X.list$X[[i]]$G[DandS]), sum(X.list$X[[i]]$G[SnotD]), sum(X.list$X[[i]]$G[DnotS]),sum(X.list$X[[i]]$G[notDS]))  
          }         
        }

    optim.out <- optim(nUS, popsize.loglik, method = "L-BFGS-B",lower=lower, upper=upper, control = list(fnscale = -1), hessian = TRUE, X=X, USdam=USdam, USsire=USsire, ped=ped, USsiredam=USsiredam)

    if(USsiredam){
      C<-diag(nbetaD+nbetaS) 
      if(nbetaD==1){
        C[1]<-solve(-1 * optim.out$hessian)
      }else{
        C[,1:nbetaD][1:nbetaD,]<-solve(-1 * optim.out$hessian)
      }
      B.start <- rep(optim.out$par,2)
    }else{
      C <- solve(-1 * optim.out$hessian)
      B.start <- optim.out$par
    }
    list(nUS=B.start, C=C)
} 
