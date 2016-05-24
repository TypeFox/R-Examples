"getsPandtP"<-function(sP, tP, PdP, GdP, X.list, nbeta, unique_id, checkP){

  if(sP$estP==TRUE & sP$estG==FALSE & length(sP$G)==0){CERVUS<-TRUE}else{CERVUS<-FALSE} 
  # use an approximation for genotyping error

  if(is.null(sP$id)){
    sP$id<-X.list$id
  }

  if(sP$estA == TRUE | sP$estG==TRUE | CERVUS==TRUE | is.null(sP$G)==FALSE){ 
    if(length(sP$A)==0){                           
      sP$A<-extractA(GdP$G)      
    }
  }

  No.E<-length(unique(GdP$categories))*length(GdP$G)*GdP$perlocus+((1-GdP$perlocus)*length(unique(GdP$categories)))

  if(is.null(sP$E1)){       
    if(is.null(GdP$categories)){    
      sP$E1<-0.005
    }else{
      sP$E1<-rep(0.005,No.E)
      if(GdP$perlocus==FALSE){
        names(sP$E1)<-unique(as.character(GdP$categories))
      }else{
        names(sP$E1)<-paste(unique(as.character(GdP$categories)), rep(names(GdP$G), each=length(unique(GdP$categories))), sep=".")
      }
    }
  }else{
    if(sP$estG==TRUE | any(sP$E1==0)){
      sP$E1[which(sP$E1==0)]<-1e-5
    }
    if(length(sP$E1)!=No.E){stop("starting values for E1 are the wrong length")}
    if(GdP$perlocus==FALSE){
      names(sP$E1)<-unique(GdP$categories)
    }else{
      names(sP$E1)<-paste(unique(GdP$categories), rep(names(GdP$G), each=length(unique(GdP$categories))), sep=".")
    }
  }
  if(GdP$marker.type=="MSC"){
     sP$estE1<-FALSE
     sP$E1<-rep(0, length(sP$E1))
  }

  if(is.null(sP$E2)){       
    if(is.null(GdP$categories)){    
      sP$E2<-0.005
    }else{
      sP$E2<-rep(0.005,No.E)
      if(GdP$perlocus==FALSE){
        names(sP$E2)<-unique(GdP$categories)
      }else{
        names(sP$E2)<-paste(unique(GdP$categories), rep(names(GdP$G), each=length(unique(GdP$categories))), sep=".")
      }
    }
  }else{
    if(sP$estG==TRUE | any(sP$E2==0)){
      sP$E2[which(sP$E2==0)]<-1e-5
    }
    if(length(sP$E2)!=No.E){stop("starting values for E1 are the wrong length")}
    if(GdP$perlocus==FALSE){
      names(sP$E2)<-unique(GdP$categories)
    }else{
      names(sP$E2)<-paste(unique(GdP$categories), rep(names(GdP$G), each=length(unique(GdP$categories))), sep=".")
    }
  }

  # if sP$G does not exist obtain from the first genotype of each individual in GdP$G

  if(is.null(sP$G)){
    sPGgiven<-FALSE
  }else{
    sPGgiven<-TRUE
  }

    if(sPGgiven == FALSE){ 
       sP$G<-lapply(GdP$G, function(x){x[-duplicated(GdP$id)==FALSE]}) 
    }else{
      for(i in 1:length(sP$A)){
        if(any((allele.names(sP$G[[i]])%in%names(sP$A[[i]]))==FALSE)){
          stop("some alleles in sP$G are not in sP$A")
        }
      }
    }
        

    if(is.null(sP$ped)==FALSE){
      if(sPGgiven){
         sP$G<-lapply(sP$G, function(x){x[order(as.numeric(match(sP$id, unique_id)))]}) 
      }
      sP$ped[,1]<-match(sP$ped[,1], unique_id)
      sP$ped[,2]<-match(sP$ped[,2], unique_id)
      sP$ped[,3]<-match(sP$ped[,3], unique_id)
      sP$ped<-sP$ped[order(as.numeric(sP$ped[,1])),]
    }

    if(sP$estUSdam==TRUE | sP$estUSsire==TRUE){
      if(sP$estUSsire=="USdam"){
        MLENus<-MLE.popsize(X.list, USdam=PdP$USdam, USsire="USdam", ped=sP$ped)
      }else{
        MLENus<-MLE.popsize(X.list, USdam=PdP$USdam, USsire=PdP$USsire, ped=sP$ped)
      }
    }

    if(length(PdP$USdam)==1 & PdP$USdam[1]==FALSE){
      nusd<-0
    }else{
      nusd<-length(unique(PdP$USdam))
    }

    if(length(PdP$USsire)==1 & PdP$USsire[1]==FALSE){
      sP$estUSsire<-FALSE
      nuss<-0
    }else{
      nuss<-length(unique(PdP$USsire))
    }

  if(sP$estUSdam==TRUE & is.null(sP$USdam)){
    sP$USdam<-MLENus$nUS[1:nusd]
  }

  if((sP$estUSsire==TRUE | sP$estUSsire=="USdam") & is.null(sP$USsire)){
    sP$USsire<-MLENus$nUS[(nusd+1):(nusd+nuss)]
  }

 if(is.null(sP$ped)){
    if(sP$estP==TRUE){  
       sP$ped<-MLE.ped(X.list, ped=NULL, USdam=PdP$USdam, nUSdam=sP$USdam, USsire=PdP$USsire, nUSsire=sP$USsire, checkP=checkP)$P   
    }else{
       sP$ped<-matrix(NA, length(sP$id), 3)
       sP$ped[,1]<-match(sP$id, unique_id)
    }
  }

  if(sP$estbeta==TRUE){
    MLEestimates<-MLE.beta(X.list, ped=sP$ped, beta=sP$beta, nUSdam=sP$USdam, nUSsire=sP$USsire)
  }
  
 
  if(sP$estbeta==TRUE){
    if(is.null(sP$beta)){
      sP$beta<-MLEestimates$beta
    }else{
      if(length(sP$beta)>1){
        sP$beta<-sP$beta[X.list$beta_map]
      }
    }
  }


  if(is.null(GdP$G)==FALSE & CERVUS==FALSE){
    if(sPGgiven==TRUE){
      if(is.null(PdP$timevar)){
        time_born=NULL
      }else{
        time_born = PdP$timevar[which(PdP$offspring==1)][match(sP$ped[,1], PdP$id)] 
      }
      if(legalG(sP$G, sP$A, sP$ped, time_born=time_born, marker.type=GdP$marker.type)$valid=="FALSE"){
        stop("sP$G does not have postive probability given possible starting pedigree")
      }
    }else{   
      if(is.null(GdP$G)==FALSE){   
        sP$G<-legalG(sP$G, sP$A, sP$ped, marker.type=GdP$marker.type)$G
      }  
    }
  }
  # get a legal configuration for sP$G 

############################################ get tuning parameters ###############################################

      if(sP$estE1){
        if(is.null(tP$E1)){  
          tP$E1<-chol(diag(No.E)*0.00003)
        }else{
          tP$E1<-chol(diag(No.E)*0.00003*tP$E1)
        }
      }

      if(sP$estE2==TRUE){
        if(is.null(tP$E2)){
          tP$E2<-chol(diag(No.E)*0.00003)
        }else{
          tP$E2<-chol(diag(No.E)*0.00003*tP$E2)
        }
      }


    if(sP$estbeta==TRUE){
      if(is.null(tP$beta)){  
        tP$beta<-rep(sqrt(10),length(sP$beta))
      }else{
        if(length(tP$beta)>1){
          tP$beta<-tP$beta[X.list$beta_map]*rep(sqrt(10),length(sP$beta))
        }else{
          tP$beta<-tP$beta*rep(sqrt(10),length(sP$beta))
        }
      }
      tP$beta<-sqrt(tP$beta%*%t(tP$beta))*MLEestimates$C
      tP$beta<-chol(tP$beta)
    }

    if(sP$estUSdam==TRUE | sP$estUSsire==TRUE){
       if(sum(diag(MLENus$C)<0)>0){
         warning("Hessian not positive-definite for MLE.popsize")
       }
    }

    if(sP$estUSdam | (nusd>0 & sP$estUSsire)){
      if(is.null(tP$USdam)){
        tP$USdam<-rep(10,nusd)
      }else{
        tP$USdam<-tP$USdam*rep(10, nusd)
      } 
      tP$USdam<-abs(tP$USdam*diag(MLENus$C)[1:nusd])
    }

     if(sP$estUSsire==TRUE | sP$estUSsire=="USdam" | (nuss>0 & sP$estUSdam)){
       if(is.null(tP$USsire)){
         tP$USsire<-rep(10, nuss)
       }else{
         tP$USsire<-tP$USsire*rep(10, nuss)
       }
       tP$USsire<-abs(tP$USsire*diag(MLENus$C)[nusd+(1:nuss)])
    }
list(sP=sP, tP=tP)
}
