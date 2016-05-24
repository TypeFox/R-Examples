getXlist<-function(PdP, GdP=NULL, A=NULL, E1=0.005, E2=0.005, mm.tol=999){

  if(is.null(GdP$id)==FALSE & is.null(PdP$id)==FALSE){
    if(FALSE%in%(GdP$id%in%PdP$id)){
      stop("genotype data exists for individuals not in PdataPed object")
    }
    if(FALSE%in%(PdP$id%in%GdP$id)){
      stop("some individuals in PdataPed object have no genotype data: replace with NA")
    }
  }

  if(is.null(PdP$id)){
    X.list<-list(id=NULL) 
    unique_id<-as.character(unique(GdP$id))
    X.list$id<-unique_id
  }else{

  null_mat<-t(as.matrix(as.numeric(NULL)))

  X.list<-list(id=NULL,beta_map=NULL, merge=c(), mergeUS=c(), X=lapply(PdP$id[which(PdP$offspring==1)], function(x){x=list(dam.id=NULL,   sire.id=NULL, mergeN=matrix(NA,2,0), XDus=null_mat, vtDus=NULL, XDs=null_mat, vtDs=NULL, XSus=null_mat, vtSus=NULL, XSs=null_mat, vtSs=NULL, XDSus=null_mat, vtDSus=NULL,XDSs=null_mat, vtDSs=NULL, G=NULL)}))   

  unique_id<-as.character(unique(PdP$id))

  X.list$id<-unique_id

  PdP$id<-match(PdP$id, unique_id)  # convert phenotypic id's to numeric

  if(length(PdP$USdam)!=1 | PdP$USdam[1]!=FALSE){
    PdP$id<-c(PdP$id, length(unique_id)+1) 
    if(is.null(PdP$sex)==FALSE){
      PdP$sex<-as.factor(c(as.character(PdP$sex), "Female")) 
    } 
    ud<-TRUE
  }else{
    ud<-FALSE
  }

  if(length(PdP$USsire)!=1 | PdP$USsire[1]!=FALSE){
    PdP$id<-c(PdP$id, length(unique_id)+ud+1)
    if(is.null(PdP$sex)==FALSE){
      PdP$sex<-as.factor(c(as.character(PdP$sex), "Male")) 
    }
    us<-TRUE
  }else{
    us<-FALSE
  }


  data_us<-matrix(NA, ud+us, length(PdP$data[1,]))
  PdP$timevar<-c(PdP$timevar, rep(NA, ud+us))
  colnames(data_us)<-colnames(PdP$data)
  PdP$data<-rbind(PdP$data, data_us)

  names(X.list$X)<-PdP$id[which(PdP$offspring==1)]

  findrest<-function(x){ # function for finding restriction variables
     if(length(grep("restrict *= *NULL" , as.character(x)))==0 & length(grep("restrict" , as.character(x)))!=0){
      int<-1
     }else{
      int<-0
     }
     int
  } 

  restrictions<-which(unlist(lapply(PdP$formula, findrest))==1)
  main_effects<-which(unlist(lapply(PdP$formula, length))==1)
  main_effects<-main_effects[main_effects%in%restrictions==FALSE]
  interactions<-which(unlist(lapply(PdP$formula, length))==2)
  interactions<-interactions[interactions%in%restrictions==FALSE]

  tmain_effects<-length(main_effects)

  if(length(interactions)>0){
    for(i in 1:length(interactions)){
      form.comb<-match(PdP$formula[[interactions[i]]], PdP$formula[main_effects[1:tmain_effects]])
      if(any(is.na(form.comb))){
        main_effects<-c(main_effects, length(PdP$formula)+1:sum(is.na(form.comb)))
        PdP$formula[length(PdP$formula)+1:sum(is.na(form.comb))]<-PdP$formula[[interactions[i]]][which(is.na(form.comb))]
      }
    }
  }
 
for(off in 1:sum(PdP$offspring==1)){

  PdP$off_record<-which(PdP$offspring==1)[off]
  PdP$keepDam<-unique(PdP$id)
  PdP$keepSire<-unique(PdP$id)
  PdP$restDam<-unique(PdP$id)
  PdP$restSire<-unique(PdP$id)

  predictors<-lapply(PdP$formula[restrictions], eval, envir=PdP)

  if(length(predictors)!=0){  
    for(i in 1:length(predictors)){
      PdP$keepDam<-PdP$keepDam[which(PdP$keepDam%in%predictors[[i]]$Dam$id==TRUE)]
      PdP$keepSire<-PdP$keepSire[which(PdP$keepSire%in%predictors[[i]]$Sire$id==TRUE)]
      PdP$restDam<-PdP$restDam[which(PdP$restDam%in%predictors[[i]]$Dam_restrict$id==TRUE)]
      PdP$restSire<-PdP$restSire[which(PdP$restSire%in%predictors[[i]]$Sire_restrict$id==TRUE)]
    }
  }else{
    if(length(PdP$sex)>0){
      PdP$keepDam<-unique(PdP$keepDam[which(PdP$sex=="Female")])
      PdP$keepSire<-unique(PdP$keepSire[which(PdP$sex=="Male")])
      PdP$restDam<-unique(PdP$restDam[which(PdP$sex=="Female")])
      PdP$restSire<-unique(PdP$restSire[which(PdP$sex=="Male")])
    }
  }

  predictors<-lapply(PdP$formula[main_effects], eval, envir=PdP)

  nvar<-rep(0, 6)         # no parameters

  if(length(predictors)!=0){  
    for(i in 1:tmain_effects){ # itterate through variables
      if(length(predictors[[i]]$Dam$X)!=0){
        nvar[1]<-nvar[1]+sum(is.na(colSums(predictors[[i]]$Dam$X)))         # starting column no. for each dam factor
        nvar[2]<-nvar[2]+sum(is.na(colSums(predictors[[i]]$Dam$X))==FALSE)  # starting column no. for each dam factor
      }
      if(length(predictors[[i]]$Sire$X)!=0){
        nvar[3]<-nvar[3]+sum(is.na(colSums(predictors[[i]]$Sire$X)))         # starting column no. for each dam factor
        nvar[4]<-nvar[4]+sum(is.na(colSums(predictors[[i]]$Sire$X))==FALSE)  # starting column no. for each dam factor
      }
      if(length(predictors[[i]]$DamSire$X)!=0){
        nvar[5]<-nvar[5]+sum(is.na(colSums(predictors[[i]]$DamSire$X)))   
        nvar[6]<-nvar[6]+sum(is.na(colSums(predictors[[i]]$DamSire$X))==FALSE)  
      }
    }
  }

  nbeta<-sum(nvar)
  X.list$X[[off]]$dam.id<-PdP$keepDam
  X.list$X[[off]]$sire.id<-PdP$keepSire
  X.list$X[[off]]$restdam.id<-PdP$restDam
  X.list$X[[off]]$restsire.id<-PdP$restSire

  ndam<-length(X.list$X[[off]]$dam.id)
  nsire<-length(X.list$X[[off]]$sire.id)

  if(nvar[1]>0){
    X.list$X[[off]]$XDus<-matrix(NA, ndam, nvar[1])
    colnames(X.list$X[[off]]$XDus)<-rep("G", nvar[1])
    X.list$X[[off]]$vtDus<-rep(NA, nvar[1])
  }
  if(nvar[2]>0){
    X.list$X[[off]]$XDs<-matrix(NA, ndam, nvar[2])
    colnames(X.list$X[[off]]$XDs)<-rep("G", nvar[2])
    X.list$X[[off]]$vtDs<-rep(NA, nvar[2])
  }
  if(nvar[3]>0){
    X.list$X[[off]]$XSus<-matrix(NA, nsire, nvar[3])
    colnames(X.list$X[[off]]$XSus)<-rep("G", nvar[3])
    X.list$X[[off]]$vtSus<-rep(NA, nvar[3])
  }
  if(nvar[4]>0){
    X.list$X[[off]]$XSs<-matrix(NA, nsire, nvar[4])
    colnames(X.list$X[[off]]$XSs)<-rep("G", nvar[4])
    X.list$X[[off]]$vtSs<-rep(NA, nvar[4])
  }
  if(nvar[5]>0){ 
    X.list$X[[off]]$XDSus<-matrix(NA, ndam*nsire, nvar[5])
    colnames(X.list$X[[off]]$XDSus)<-rep("G",nvar[5])
    X.list$X[[off]]$vtDSus<-rep(NA, nvar[5])
  }
  if(nvar[6]>0){ 
    X.list$X[[off]]$XDSs<-matrix(NA, ndam*nsire, nvar[6])
    colnames(X.list$X[[off]]$XDSs)<-rep("G",nvar[6])
    X.list$X[[off]]$vtDSs<-rep(NA, nvar[6])
  }
  # sets up empty design matrix ncolumns = npredictors+1 for genetic likelihoods

  ##########################################################################################################
  ###################################### main effects ######################################################
  ##########################################################################################################

  if(tmain_effects!=0){  

    nvar_tmp<-rep(0,6) 
 
    for(i in 1:tmain_effects){ # itterates through the variables
      # Dam variables 
      if(length(predictors[[i]]$Dam$X)!=0){ 
        if(is.na(sum(predictors[[i]]$Dam$X))==TRUE){
          for(c in 1:ncol(predictors[[i]]$Dam$X)){
            nvar_tmp[1]<-nvar_tmp[1]+1       
            X.list$X[[off]]$vtDus[nvar_tmp[1]]<-predictors[[i]]$Dam$var_type
            X.list$X[[off]]$XDus[,nvar_tmp[1]]<-predictors[[i]]$Dam$X[,c] 
            colnames(X.list$X[[off]]$XDus)[nvar_tmp[1]]<-predictors[[i]]$Dam$var_name[c]
            if(any(is.na(X.list$X[[off]]$XDus[,nvar_tmp[1]][-ndam]))){stop("Missing covariate data")}
            if(predictors[[i]]$Dam$merge==TRUE){
              if(off==1){
                X.list$merge<-c(X.list$merge, nvar_tmp[1])
                X.list$mergeUS<-c(X.list$mergeUS, 0)
              }
              X.list$X[[off]]$mergeN<-cbind(X.list$X[[off]]$mergeN, c(sum(predictors[[i]]$Dam$X[,c]==1, na.rm=T), sum(predictors[[i]]$Dam$X[,c]==0, na.rm=T))) 
            }
          }
        }else{
          for(c in 1:ncol(predictors[[i]]$Dam$X)){
            nvar_tmp[2]<-nvar_tmp[2]+1       
            X.list$X[[off]]$vtDs[nvar_tmp[2]]<-predictors[[i]]$Dam$var_type
            X.list$X[[off]]$XDs[,nvar_tmp[2]]<-predictors[[i]]$Dam$X[,c] 
            colnames(X.list$X[[off]]$XDs)[nvar_tmp[2]]<-predictors[[i]]$Dam$var_name[c]
            if(any(is.na(X.list$X[[off]]$XDs[,nvar_tmp[2]]))){stop("Missing covariate data")}
            if(predictors[[i]]$Dam$merge==TRUE){
              if(off==1){
                X.list$merge<-c(X.list$merge, nvar[1]+nvar_tmp[2])
                X.list$mergeUS<-c(X.list$mergeUS, ud*((predictors[[i]]$Dam$X[,c][nrow(predictors[[i]]$Dam$X)]==0)+1))
              }
              X.list$X[[off]]$mergeN<-cbind(X.list$X[[off]]$mergeN, c(sum(predictors[[i]]$Dam$X[,c]==1), sum(predictors[[i]]$Dam$X[,c]==0))) 
            }
          }
        }
      }


      #Sire variables 
      if(length(predictors[[i]]$Sire$X)!=0){ 
        if(is.na(sum(predictors[[i]]$Sire$X))==TRUE){
          for(c in 1:ncol(predictors[[i]]$Sire$X)){
            nvar_tmp[3]<-nvar_tmp[3]+1
            X.list$X[[off]]$vtSus[nvar_tmp[3]]<-predictors[[i]]$Sire$var_type
            X.list$X[[off]]$XSus[,nvar_tmp[3]]<-predictors[[i]]$Sire$X[,c] 
            colnames(X.list$X[[off]]$XSus)[nvar_tmp[3]]<-predictors[[i]]$Sire$var_name[c]
            if(any(is.na(X.list$X[[off]]$XSus[,nvar_tmp[3]][-nsire]))){stop("Missing covariate data")}
            if(predictors[[i]]$Sire$merge==TRUE){
              if(off==1){
                X.list$merge<-c(X.list$merge, sum(nvar[1:2])+nvar_tmp[3])
                X.list$mergeUS<-c(X.list$mergeUS, 0)
              }
              X.list$X[[off]]$mergeN<-cbind(X.list$X[[off]]$mergeN, c(sum(predictors[[i]]$Sire$X[,c]==1, na.rm=T), sum(predictors[[i]]$Sire$X[,c]==0, na.rm=T))) 
            }
          }
        }else{
          for(c in 1:ncol(predictors[[i]]$Sire$X)){
            nvar_tmp[4]<-nvar_tmp[4]+1
            X.list$X[[off]]$vtSs[nvar_tmp[4]]<-predictors[[i]]$Sire$var_type
            X.list$X[[off]]$XSs[,nvar_tmp[4]]<-predictors[[i]]$Sire$X[,c] 
            colnames(X.list$X[[off]]$XSs)[nvar_tmp[4]]<-predictors[[i]]$Sire$var_name[c]
            if(any(is.na(X.list$X[[off]]$XSs[,nvar_tmp[4]]))){stop("Missing covariate data")}
            if(predictors[[i]]$Sire$merge==TRUE){
              if(off==1){
                X.list$merge<-c(X.list$merge, sum(nvar[1:3])+nvar_tmp[4])
                X.list$mergeUS<-c(X.list$mergeUS, us*((predictors[[i]]$Sire$X[,c][nrow(predictors[[i]]$Sire$X)]==0)+1))
              }
              X.list$X[[off]]$mergeN<-cbind(X.list$X[[off]]$mergeN, c(sum(predictors[[i]]$Sire$X[,c]==1, na.rm=T), sum(predictors[[i]]$Sire$X[,c]==0, na.rm=T))) 
            }
          }
        }
      }

    #Dam/Sire variables 
      if(length(predictors[[i]]$DamSire$X)!=0){
        if(is.na(sum(predictors[[i]]$DamSire$X))==TRUE){
          for(c in 1:ncol(predictors[[i]]$DamSire$X)){
            nvar_tmp[5]<-nvar_tmp[5]+1
            X.list$X[[off]]$vtDSus[nvar_tmp[5]]<-predictors[[i]]$DamSire$var_type
            X.list$X[[off]]$XDSus[,nvar_tmp[5]]<-predictors[[i]]$DamSire$X[,c]
            colnames(X.list$X[[off]]$XDSus)[nvar_tmp[5]]<-predictors[[i]]$DamSire$var_name[c]
            if(us==TRUE){rem.var<-seq(nsire,ndam*nsire, nsire)}
            if(ud==TRUE){rem.var<-((((ndam-1)*nsire)+1):(ndam*nsire))}
            if(us==TRUE & ud==TRUE){rem.var<-c(seq(nsire,ndam*nsire, nsire), (((ndam-1)*nsire)+1):c((ndam*nsire)-1))}
            if(any(is.na(X.list$X[[off]]$XDSus[,nvar_tmp[5]][-rem.var]))){stop("Missing covariate data")}
          }
        }else{
          for(c in 1:ncol(predictors[[i]]$DamSire$X)){
            nvar_tmp[6]<-nvar_tmp[6]+1
            X.list$X[[off]]$vtDSs[nvar_tmp[6]]<-predictors[[i]]$DamSire$var_type
            X.list$X[[off]]$XDSs[,nvar_tmp[6]]<-predictors[[i]]$DamSire$X[,c]
            colnames(X.list$X[[off]]$XDSs)[nvar_tmp[6]]<-predictors[[i]]$DamSire$var_name[c]     
            if(any(is.na(X.list$X[[off]]$XDSs[,nvar_tmp[6]]))){stop("Missing covariate data")}        
           }
        }
      }
    }
  }
  
  ###################################################################################################################
  ##################################  interactions ##################################################################
  ###################################################################################################################

  par_type=c()
  if(length(interactions)>0){
    
    for(i in 1:length(interactions)){ 

      form.comb<-match(PdP$formula[[interactions[i]]], PdP$formula[main_effects])
      t1<-predictors[[form.comb[1]]]
      t2<-predictors[[form.comb[2]]]

      if(off==1){

        if(i==1){
          par_type=c()
          dam.dam=rep(FALSE, length(interactions))
          sire.sire=rep(FALSE, length(interactions))
          dam.sire=rep(FALSE, length(interactions))
          damsire.damsire=rep(FALSE, length(interactions))
          effect=rep(0, length(interactions))
        }
        if(is.null(t1$Dam$X)==FALSE & is.null(t1$Sire$X)==FALSE){
          if(is.null(t2$Dam$X)==FALSE & is.null(t2$Sire$X)==FALSE){
            dam.dam[i]=TRUE
            sire.sire[i]=TRUE
          }else{
            stop("interactions between a genderless variable and a sex-specific variable not possible")
          }
        }     

        if(is.null(t1$Dam$X)==FALSE & is.null(t2$Dam$X)==FALSE){
          dam.dam[i]=TRUE
          if(TRUE%in%(is.na(t1$Dam$X)) | TRUE%in%(is.na(t2$Dam$X))){
            effect[i]<-1
            par_type<-c(par_type, rep(effect[i], ncol(t1$Dam$X)*ncol(t2$Dam$X)))
          }else{
            effect[i]<-2
            par_type<-c(par_type, rep(effect[i], ncol(t1$Dam$X)*ncol(t2$Dam$X)))
          }
        }
        
        if(is.null(t1$Sire$X)==FALSE & is.null(t2$Sire$X)==FALSE){
          sire.sire[i]=TRUE
          if(TRUE%in%(is.na(t1$Sire$X)) | TRUE%in%(is.na(t2$Sire$X))){
            effect[i]<-3
            par_type<-c(par_type, rep(effect[i], ncol(t1$Sire$X)*ncol(t2$Sire$X)))
          }else{
            effect[i]<-4
            par_type<-c(par_type, rep(effect[i], ncol(t1$Sire$X)*ncol(t2$Sire$X)))
          }
        }

        if(is.null(t1$Dam$X)==FALSE & is.null(t2$Sire$X)==FALSE){
          if(is.null(t2$Dam$X) & is.null(t1$Sire$X)){
            dam.sire[i]=TRUE
            if(TRUE%in%(is.na(t1$Dam$X)) | TRUE%in%(is.na(t2$Sire$X))){
              effect[i]<-5
              par_type<-c(par_type, rep(effect[i], ncol(t1$Dam$X)*ncol(t2$Sire$X)))
            }else{
              effect[i]<-6
              par_type<-c(par_type, rep(effect[i], ncol(t1$Dam$X)*ncol(t2$Sire$X)))
            }
          }
        }

        if(is.null(t2$Dam$X)==FALSE & is.null(t1$Sire$X)==FALSE){
          if(is.null(t1$Dam$X) & is.null(t2$Sire$X)){
            dam.sire[i]=TRUE
            if(TRUE%in%(is.na(t1$Sire$X)) | TRUE%in%(is.na(t2$Dam$X))){
              effect[i]<-5
              par_type<-c(par_type, rep(effect[i], ncol(t1$Sire$X)*ncol(t2$Dam$X)))
            }else{
              effect[i]<-6
              par_type<-c(par_type, rep(effect[i], ncol(t1$Sire$X)*ncol(t2$Dam$X)))
            }
          }
        }

        if(is.null(t1$DamSire$X)==FALSE & is.null(t2$DamSire$X)==FALSE){
          damsire.damsire[i]=TRUE
            if(TRUE%in%(is.na(t1$DamSire$X)) | TRUE%in%(is.na(t2$DamSire$X))){
              effect[i]<-5
              par_type<-c(par_type, rep(effect[i], ncol(t1$DamSire$X)*ncol(t2$DamSire$X)))
            }else{
              effect[i]<-6
              par_type<-c(par_type, rep(effect[i], ncol(t1$DamSire$X)*ncol(t2$DamSire$X)))
            }
        }
      }

      col<-0 
      if(dam.dam[i]){
        int.tmp<-matrix(NA,nrow(t1$Dam$X), ncol(t1$Dam$X)*ncol(t2$Dam$X))
        colnames(int.tmp)<-rep("G", ncol(int.tmp))
        for(v1 in 1:ncol(t1$Dam$X)){
          for(v2 in 1:ncol(t2$Dam$X)){
            col<-col+1
            int.tmp[,col]<-t1$Dam$X[,v1]*t2$Dam$X[,v2]
            colnames(int.tmp)[col]<-paste(t1$Dam$var_name[v1], t2$Dam$var_name[v2], sep=".")
         }
        }
        if(effect[i]==1){
          for(c in 1:ncol(int.tmp)){
            nvar[1]<-nvar[1]+1       
            X.list$X[[off]]$XDus<-as.matrix(cbind(X.list$X[[off]]$XDus, int.tmp[,c]))
            if(t1$Dam$var_type == "factor" & t2$Dam$var_type == "factor"){
              X.list$X[[off]]$vtDus<-c(X.list$X[[off]]$vtDus, "factor")
            }else{
              X.list$X[[off]]$vtDus<-c(X.list$X[[off]]$vtDus, "numeric")
            }
            colnames(X.list$X[[off]]$XDus)[nvar[1]]<-colnames(int.tmp)[c]
          }
        }else{
          for(c in 1:ncol(int.tmp)){
            nvar[2]<-nvar[2]+1       
            X.list$X[[off]]$XDs<-as.matrix(cbind(X.list$X[[off]]$XDs, int.tmp[,c]))
            if(t1$Dam$var_type == "factor" & t2$Dam$var_type == "factor"){
              X.list$X[[off]]$vtDs<-c(X.list$X[[off]]$vtDs, "factor")
            }else{
              X.list$X[[off]]$vtDs<-c(X.list$X[[off]]$vtDs, "numeric")
            }
            colnames(X.list$X[[off]]$XDs)[nvar[2]]<-colnames(int.tmp)[c]
          }
        }
      }

      col<-0 
      if(sire.sire[i]){
        int.tmp<-matrix(NA,nrow(t1$Sire$X), ncol(t1$Sire$X)*ncol(t2$Sire$X))
        colnames(int.tmp)<-rep("G", ncol(int.tmp))
        for(v1 in 1:ncol(t1$Sire$X)){
          for(v2 in 1:ncol(t2$Sire$X)){
            col<-col+1
            int.tmp[,col]<-t1$Sire$X[,v1]*t2$Sire$X[,v2]
            colnames(int.tmp)[col]<-paste(t1$Sire$var_name[v1], t2$Sire$var_name[v2], sep=".")
          }
        }
        if(effect[i]==3){
          for(c in 1:ncol(int.tmp)){
            nvar[3]<-nvar[3]+1       
            X.list$X[[off]]$XSus<-as.matrix(cbind(X.list$X[[off]]$XSus, int.tmp[,c]))
            if(t1$Sire$var_type == "factor" & t2$Sire$var_type == "factor"){
              X.list$X[[off]]$vtSus<-c(X.list$X[[off]]$vtSus, "factor")
            }else{
              X.list$X[[off]]$vtSus<-c(X.list$X[[off]]$vtSus, "numeric")
            }
            colnames(X.list$X[[off]]$XSus)[nvar[3]]<-colnames(int.tmp)[c]
          }
        }else{
          for(c in 1:ncol(int.tmp)){
            nvar[4]<-nvar[4]+1       
            X.list$X[[off]]$XSs<-as.matrix(cbind(X.list$X[[off]]$XSs, int.tmp[,c]))
            if(t1$Sire$var_type == "factor" & t2$Sire$var_type == "factor"){
              X.list$X[[off]]$vtSs<-c(X.list$X[[off]]$vtSs, "factor")
            }else{
              X.list$X[[off]]$vtSs<-c(X.list$X[[off]]$vtSs, "numeric")
            }
            colnames(X.list$X[[off]]$XSs)[nvar[4]]<-colnames(int.tmp)[c]
          }
        }
      }

      col<-0 
      if(dam.sire[i] | damsire.damsire[i]){
        if(dam.sire[i]){
          int.tmp<-matrix(NA,nrow(t1$Dam$X), ncol(t1$Dam$X)*ncol(t2$Sire$X))
          colnames(int.tmp)<-rep("G", ncol(int.tmp))
           for(v1 in 1:ncol(t1$Dam$X)){
            for(v2 in 1:ncol(t2$Sire$X)){
              col<-col+1
              nsires<-length(X.list$X[[off]]$sire.id)
              ndams<-length(X.list$X[[off]]$dam.id)
              int.tmp[,col]<-rep(t1$Dam$X[,v1], each=nsires)*rep(t2$Sire$X[,v2], ndams)
              colnames(int.tmp)[col]<-paste(t1$Sire$var_name[v1], t2$Sire$var_name[v2], sep=".")
            }
          }
        }else{
          int.tmp<-matrix(NA,nrow(t1$DamSire$X), ncol(t1$DamSire$X)*ncol(t2$DamSire$X))
          colnames(int.tmp)<-rep("G", ncol(int.tmp))
          for(v1 in 1:ncol(t1$DamSire$X)){
            for(v2 in 1:ncol(t2$DamSire$X)){
              col<-col+1
              int.tmp[,col]<-t1$DamSire$X[,v1]*t2$DamSire$X[,v2]
              colnames(int.tmp)[col]<-paste(t1$DamSire$var_name[v1], t2$DamSire$var_name[v2], sep=".")
            }
          }
        }
        if(effect[i]==5){
          if(ncol(X.list$X[[off]]$XDSus)==0){
            X.list$X[[off]]$XDSus<-matrix(0, length(X.list$X[[off]]$sire.id)*length(X.list$X[[off]]$dam.id), 0)
          }
          for(c in 1:ncol(int.tmp)){
            nvar[5]<-nvar[5]+1       
            X.list$X[[off]]$XDSus<-as.matrix(cbind(X.list$X[[off]]$XDSus, int.tmp[,c]))
            if(dam.sire[i]){
              if(t1$Dam$var_type == "factor" & t2$Sire$var_type == "factor"){
                X.list$X[[off]]$vtDSus<-c(X.list$X[[off]]$vtDSus, "factor")
              }else{
                X.list$X[[off]]$vtDSus<-c(X.list$X[[off]]$vtDSus, "numeric")
              }
            }else{
              if(t1$DamSire$var_type == "factor" & t2$DamSire$var_type == "factor"){
                X.list$X[[off]]$vtDSus<-c(X.list$X[[off]]$vtDSus, "factor")
              }else{
                X.list$X[[off]]$vtDSus<-c(X.list$X[[off]]$vtDSus, "numeric")
              }
            }
            colnames(X.list$X[[off]]$XDSus)[nvar[5]]<-colnames(int.tmp)[c]
          }
        }else{
          if(ncol(X.list$X[[off]]$XDSs)==0){
            X.list$X[[off]]$XDSs<-matrix(0, length(X.list$X[[off]]$sire.id)*length(X.list$X[[off]]$dam.id), 0)
          }
          for(c in 1:ncol(int.tmp)){
            nvar[6]<-nvar[6]+1       
            X.list$X[[off]]$XDSs<-as.matrix(cbind(X.list$X[[off]]$XDSs, int.tmp[,c]))
            if(dam.sire[i]){
              if(t1$Dam$var_type == "factor" & t2$Sire$var_type == "factor"){
                X.list$X[[off]]$vtDSs<-c(X.list$X[[off]]$vtDSs, "factor")
              }else{
                X.list$X[[off]]$vtDSs<-c(X.list$X[[off]]$vtDSs, "numeric")
              }
            }else{
              if(t1$DamSire$var_type == "factor" & t2$DamSire$var_type == "factor"){
                X.list$X[[off]]$vtDSs<-c(X.list$X[[off]]$vtDSs, "factor")          
              }else{
                X.list$X[[off]]$vtDSs<-c(X.list$X[[off]]$vtDSs, "numeric")
              }
            }
          }
          colnames(X.list$X[[off]]$XDSs)[nvar[6]]<-colnames(int.tmp)[c]
        }
      }
    }
  }
}

  if(sum(nvar)>0){
    beta_map<-1:sum(nvar)
    if(sum(nvar[3:4])>0){
      Dlinked<-c(grep("linked", colnames(X.list$X[[1]]$XDus)), grep("linked", colnames(X.list$X[[1]]$XDs))+nvar[1])
      Dlinked_names<-c(colnames(X.list$X[[1]]$XDus), colnames(X.list$X[[1]]$XDs))[Dlinked]
      Slinked<-match(c(colnames(X.list$X[[1]]$XSus), colnames(X.list$X[[1]]$XSs)), Dlinked_names)
      Slinked[which(is.na(Slinked)==FALSE)]<-Dlinked
      Slinked[which(is.na(Slinked)==TRUE)]<-sum(nvar[1:2])+c(1:sum(is.na(Slinked)))
      beta_map[sum(nvar[1:2])+(1:sum(nvar[3:4]))]<-Slinked
    }
    if(sum(nvar[5:6])>0 & sum(nvar[1:4])>0){
      beta_map[sum(nvar[1:4])+(1:sum(nvar[5:6]))]<-c(max(beta_map[1:sum(nvar[1:4])])+(1:sum(nvar[5:6])))
    }
  }else{
    beta_map<--999
  }

  X.list$beta_map<-beta_map
# contrast with base parents

    for(off in 1:sum(PdP$offspring==1)){ 

      if(is.null(X.list$merge)==FALSE){
        for(m in 1:length(X.list$merge)){
          X.list$X[[off]]$mergeN[,m][X.list$mergeUS[m]]<-X.list$X[[off]]$mergeN[,m][X.list$mergeUS[m]]-1  
          # need to take 1 off the mergeN class as it is actually unsampled
          n1<-X.list$X[[off]]$mergeN[,m][1]+(X.list$mergeUS[m]==1)
          n2<-X.list$X[[off]]$mergeN[,m][2]+(X.list$mergeUS[m]==2)
          if(n1==0 | n2==0){
            X.list$X[[off]]$mergeN[,m]<-1 
          }
          # if all individuals (sampled and unsampled are in 1 class numerical problems occur) 
          # however mergeN can be safley replaced with what ever since they don't contribute 
          # to the likelihood or pedigree estimation as all individuals are monomorphic!        
        }
      }

      if(nvar[1]>0){
        nrowX=dim(X.list$X[[off]]$XDus)[1]
        ncolX=dim(X.list$X[[off]]$XDus)[2]
        base<-X.list$X[[off]]$XDus[1,]
        X.list$X[[off]]$XDus<-X.list$X[[off]]$XDus-matrix(rep(base,each=nrowX), nrowX, ncolX)  
        col2scale<-which(X.list$X[[off]]$vtDus=="numeric")
        if(length(col2scale)>0){
          center.val<-colMeans(as.matrix(X.list$X[[off]]$XDus[,col2scale]), na.rm=T)
          X.list$X[[off]]$XDus[,col2scale]<-scale(X.list$X[[off]]$XDus[,col2scale], center=center.val, scale=FALSE)
        }
      }
      if(nvar[2]>0){
        nrowX=dim(X.list$X[[off]]$XDs)[1]
        ncolX=dim(X.list$X[[off]]$XDs)[2]
        base<-X.list$X[[off]]$XDs[1,]
        X.list$X[[off]]$XDs<-X.list$X[[off]]$XDs-matrix(rep(base,each=nrowX), nrowX, ncolX) 
        col2scale<-which(X.list$X[[off]]$vtDs=="numeric")
        if(length(col2scale)>0){
          center.val<-colMeans(as.matrix(X.list$X[[off]]$XDs[,col2scale]), na.rm=T)
          X.list$X[[off]]$XDs[,col2scale]<-scale(X.list$X[[off]]$XDs[,col2scale], center=center.val, scale=FALSE)
        }
      }
      if(nvar[3]>0){
        nrowX=dim(X.list$X[[off]]$XSus)[1]
        ncolX=dim(X.list$X[[off]]$XSus)[2]
        base<-X.list$X[[off]]$XSus[1,]
        X.list$X[[off]]$XSus<-X.list$X[[off]]$XSus-matrix(rep(base,each=nrowX), nrowX, ncolX) 
        col2scale<-which(X.list$X[[off]]$vtSus=="numeric")
        if(length(col2scale)>0){
          center.val<-colMeans(as.matrix(X.list$X[[off]]$XSus[,col2scale]), na.rm=T)
          X.list$X[[off]]$XSus[,col2scale]<-scale(X.list$X[[off]]$XSus[,col2scale], center=center.val, scale=FALSE)
        } 
      }
      if(nvar[4]>0){
        nrowX=dim(X.list$X[[off]]$XSs)[1]
        ncolX=dim(X.list$X[[off]]$XSs)[2]
        base<-X.list$X[[off]]$XSs[1,]
        X.list$X[[off]]$XSs<-X.list$X[[off]]$XSs-matrix(rep(base,each=nrowX), nrowX, ncolX)
        col2scale<-which(X.list$X[[off]]$vtSs=="numeric")
        if(length(col2scale)>0){
          center.val<-colMeans(as.matrix(X.list$X[[off]]$XSs[,col2scale]), na.rm=T)
          X.list$X[[off]]$XSs[,col2scale]<-scale(X.list$X[[off]]$XSs[,col2scale], center=center.val, scale=FALSE)
        }   
      }
      if(nvar[5]>0){
        nrowX=dim(X.list$X[[off]]$XDSus)[1]
        ncolX=dim(X.list$X[[off]]$XDSus)[2]
        base<-X.list$X[[off]]$XDSus[1,]
        X.list$X[[off]]$XDSus<-X.list$X[[off]]$XDSus-matrix(rep(base,each=nrowX), nrowX, ncolX)  
        col2scale<-which(X.list$X[[off]]$vtDSus=="numeric")
        if(length(col2scale)>0){
          center.val<-colMeans(as.matrix(X.list$X[[off]]$XDSus[,col2scale]), na.rm=T)
          X.list$X[[off]]$XDSus[,col2scale]<-scale(X.list$X[[off]]$XDSus[,col2scale], center=center.val, scale=FALSE)
        } 
      }
      if(nvar[6]>0){
        nrowX=dim(X.list$X[[off]]$XDSs)[1]
        ncolX=dim(X.list$X[[off]]$XDSs)[2]
        base<-X.list$X[[off]]$XDSs[1,]
        X.list$X[[off]]$XDSs<-X.list$X[[off]]$XDSs-matrix(rep(base,each=nrowX), nrowX, ncolX) 
        col2scale<-which(X.list$X[[off]]$vtDSs=="numeric")
        if(length(col2scale)>0){
          center.val<-colMeans(as.matrix(X.list$X[[off]]$XDSs[,col2scale]), na.rm=T)
          X.list$X[[off]]$XDSs[,col2scale]<-scale(X.list$X[[off]]$XDSs[,col2scale], center=center.val, scale=FALSE)
        }  
      }
    }    

    if(is.null(GdP$G)==FALSE){
      if(is.null(A)==TRUE){
        A<-extractA(GdP$G)
      }else{
        for(i in 1:length(GdP$G)){
          A[[i]]<-A[[i]][order(A[[i]], decreasing=T)]
          GdP$G[[i]]<-genotype(GdP$G[[i]], alleles=names(A[[i]]), reorder="no")
        }
      }
      Gid<-GdP$id[-duplicated(GdP$id)==FALSE]
      G<-lapply(GdP$G, function(x){x[-duplicated(GdP$id)==FALSE]})
      grouped_by_id<-order(match(Gid, unique_id))        
      G<-lapply(G, function(x){x[grouped_by_id]}) 
      Gid<-grouped_by_id  

      X.list<-mismatches(X.list, G=G, mm.tol=mm.tol)
      if(is.null(E1)==TRUE){
        E1<-0.005
      }
      if(is.null(E2)==TRUE){
        E2<-0.005
      }

      X.list<-fillX.G(X.list, A=A, G=G, E1=E1, E2=E2, marker.type=GdP$marker.type)
      X.list<-reordXlist(X.list, marker.type=GdP$marker.type)
    }

    npdam<-unlist(lapply(X.list$X, function(x){length(x$restdam.id)}))
    npsire<-unlist(lapply(X.list$X, function(x){length(x$restsire.id)}))

    if(any(npdam==0)){ stop(paste("Indiviudals", paste(X.list$id[as.numeric(names(X.list$X)[which(npdam==0)])], collapse=" "), "have no possible dams"))}
    if(any(npsire==0)){stop(paste("Individuals", paste(X.list$id[as.numeric(names(X.list$X)[which(npsire==0)])], collapse=" "), "have no possible sires"))}

  }
X.list
}
