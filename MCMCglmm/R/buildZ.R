buildZ<-function(x, data, nginverse=NULL){

  vtype="idh"   # form of covariances for the random effets of a random term classified by some factor 
  rtype="iid"   # form of covariances between random effets of random terms

  if(length(grep("^us\\(", x))>0){
    vtype<-"us"
  }
  if(length(grep("^corh\\(", x))>0){
    vtype<-"corh"
    stop("sorry - corh not yet implemented")
  }
  if(length(grep("^cor\\(", x))>0){
    vtype<-"cor"
    stop("sorry - cor not yet implemented. Note that cor as used in versions <2.18 has now been replaced by corg")
  }
  if(length(grep("^corg\\(", x))>0){
    vtype<-"corg"
  }
  if(length(grep("^corgh\\(", x))>0){
    vtype<-"corgh"
  }
  if(length(grep("^cors\\(", x))>0){
    vtype<-"cors"
  }
  if(length(grep("^idv\\(", x))>0){
    vtype<-"idv"
  }
  if(length(grep("^ante.*\\(", x))>0){
    vtype<-gsub("\\(.*", "", x)
  }
  if(length(grep("^sub\\(", x))>0){
    vtype<-"sub"
  }
  if(length(grep("(^|:)mm\\(", x))>0){
    rtype<-"mm"
  }
  if(length(grep("(^|:)str\\(", x))>0){
    rtype<-"str"
  }

  fformula<-gsub("^(us|corg|corgh|corh|cor|cors|idh|idv|antec?[0-9]*v?|sub)\\(", "", x)

  if(grepl("^str\\(|^mm\\(", fformula)){
    fformula<-paste("):", fformula, sep="")
  }

  openB<-gregexpr("\\(", fformula)[[1]]
  closeB<-gregexpr("\\)", fformula)[[1]]

  if(openB[1]!=-1){
    while(openB[1]<closeB[1]){
      openB<-openB[-1]
      closeB<-closeB[-1]
      if(length(closeB)==0 | length(openB)==0){break}
    }
  }

  rterms<-substr(fformula,closeB[1]+2,nchar(fformula))
  rterms<-gsub("(^|:)(str|mm|)\\(", "", rterms)
  rterms<-gsub("\\)$", "", rterms)

  if(rtype!="iid" & any(grepl(":", rterms))){
    stop("interactions not permitted in str and mm structures")
  }
  if(rtype=="mm"){
   pm.pos<-gregexpr("\\+|\\-", rterms)[[1]]
   mm.sign<-c(1,sapply(pm.pos, function(x){if(substr(rterms, x,x)=="+"){1}else{-1}}))
  }

  rterms<-strsplit(rterms, ":| \\+ | \\- ")[[1]]

  if(is.null(nginverse)==FALSE){
    Aterm<-na.omit(match(rterms,nginverse))[1]
    if(is.na(Aterm)){
      Aterm<-0
    }
  }else{
    Aterm<-0
  }

  fformula<-substr(fformula,1,closeB-1)
  fformula<-gsub(" ", "", fformula)

  idv.vnames<-paste(fformula, collapse="")

  if(fformula!=""){fformula.exists=TRUE}else{fformula.exists=FALSE}

  if(fformula!=""){fformula<-c("+", fformula)}

  fformula<-as.formula(paste("~", "-1",  paste(fformula, collapse=""), sep=""))

  orig.na.action<-options("na.action")[[1]]
  options("na.action"="na.pass")

  X<-model.matrix(fformula, data)

  options("na.action"=orig.na.action)

  X<-as(X, "sparseMatrix")

  if(any(diff(X@p)==0)){
    X<-X[,-which(diff(X@p)==0)]
    # sometimes X has null colmumns (e.g. for at.level(a,1):b, where b is factor and some levels are not present when a=levels(a)[1])
  }

  nfl<-ncol(X)

  if(nfl==0 & fformula.exists){stop("variance function formula for some random term defines predictors that are all zero")}

  if(nfl==0){
    trait.ordering<-data$trait[1]
  }else{
    trait.ordering<-as.numeric(data$trait[X@i[X@p[1:ncol(X)]+1]+1])
  }

  ZZ<-list()

    for(k in 1:max(1,(1-(rtype=="iid"))*length(rterms))){  # iterate through Rterms or evalute once if iid

      if(rtype=="iid"){
        select.terms<-1:length(rterms)
      }else{
        select.terms<-k
      } 

      if(length(rterms)==0){
        rfactor<-rep(as.factor(1), dim(data)[1]) 	
      }else{
        if(any(rterms%in%colnames(data)==FALSE)){stop(paste("object", paste(rterms, collapse=" and "), "not found"))}
        rfactor<-interaction(data[,rterms[select.terms]], drop=(Aterm==0 & rtype=="iid"))
        # random effects are dropped if they are not animal or part of a mm or str structure 
      }

      nrl<-nlevels(rfactor) 
      nd<-length(rfactor)
              
      data_pos<-as.numeric(rfactor)
      ZZ[[k]]<-Matrix(0,nd,nrl)
      ZZ[[k]][,1][2]<-1              # force it out of vbeing upper triangle!!!!
      ZZ[[k]]@p<-as.integer(c(0,cumsum(table(rfactor))))    
      cnt<-0
      for(j in 1:nrl){
        hit<-which(data_pos==j)
        hit<-hit-nd*(ceiling(hit/nd)-1)
        if(length(hit)>0){
          ZZ[[k]]@i[cnt+1:length(hit)]<-as.integer(hit-1)
          cnt<-cnt+length(hit)
        }
      }
      ZZ[[k]]@x<-rep(1,length(ZZ[[k]]@i))

      colnames(ZZ[[k]])<-paste(paste(rterms[select.terms], collapse=":"),levels(rfactor), sep=".")
      if(any(data$MCMC_dummy==0 & is.na(rfactor))){warning("missing values in random predictors")}

    }

    if(nfl==0){

      nfl<-1  

      for(k in 1:max(1,(1-(rtype=="iid"))*length(rterms))){  # iterate through Rtrems or evalute once if iid

        missing<-which(diff(ZZ[[k]]@p)==0)  # non-represented random effects (or zero covariates for RR)
 
        if(length(missing)>0){        
          if(Aterm==0 & rtype=="iid"){                       # if not associated with ginv drop terms     
            ZZ[[k]]<-ZZ[[k]][,-missing,drop=FALSE]
            nrl[i]<-ncol(ZZ[k])              
          }                           
        }

        if(k>1){
          if(is.null(levels(data[,rterms[select.terms]])) | any(levels(data[,rterms[select.terms]])!=levels(data[,rterms[1]]))){
            stop("terms involved in mm/str structures must have identical levels")
          }
          if(rtype=="mm"){
            colnames(ZZ[[k]])<-colnames(ZZ[[1]])
            if(mm.sign[k]==1){
              ZZ[[1]]<-ZZ[[1]]+ZZ[[k]]
            }else{
              ZZ[[1]]<-ZZ[[1]]-ZZ[[k]]
            }
          }
          if(rtype=="str"){
            ZZ[[1]]<-cBind(ZZ[[1]],ZZ[[k]])
          }
        }
      }
      vnames<-paste(rterms, collapse=":")
      if(rtype=="mm"){
         vnames<-paste(rterms, collapse="+")
      }
      if(rtype=="str"){
        vnames<-apply(expand.grid(rterms,rterms), 1, paste, collapse=".")
        nfl<-nfl*length(rterms)
        vtype<-rep("us", length(vtype))
      }

      return(list(Z=ZZ[[1]], nfl=nfl, nrl=nrl, Aterm=Aterm, vtype=vtype, vnames=vnames, trait.ordering=trait.ordering))

    }else{

      if(vtype=="idh"){
        nrl<-rep(nrl, nfl)
      }
      
      for(i in 1:nfl){
    
        Xtmp<-Matrix(0,nd,nd)
        Xcol<-X[,i,drop=FALSE]
        Xtmp@p<-as.integer(rep(0:c(length(Xcol@i)-1),diff(c(0,Xcol@i+1))))
        Xtmp@p[(length(Xtmp@p)+1):(nd+1)]<-length(Xcol@i)
        Xtmp@i<-Xcol@i
        Xtmp@x<-Xcol@x

        Z<-list()

        for(k in 1:max(1,(1-(rtype=="iid"))*length(rterms))){  # iterate through Rtrems or evalute once if iid

          if(rtype=="iid"){
            select.terms<-1:length(rterms)
          }

          Z[[k]]<-Xtmp%*%ZZ[[k]]
          colnames(Z[[k]])<-paste(colnames(X)[i], colnames(Z[[k]]),sep=".") 

          missing<-which(diff(Z[[k]]@p)==0)

          #################################################################################################
          # It woule be more efficient to allow complete sparse columns in Z - but need to change C code! #
          #################################################################################################

          if(length(missing)>0){
            if(vtype=="idh" | vtype=="idv"){
              if(Aterm==0 & rtype=="iid"){                     # if univariate G-structure and no Ginv term drop levels
                Z[[k]]<-Z[[k]][,-missing,drop=FALSE]
                nrl[i]<-ncol(Z[[k]])
              }
            }
          }

          if(k>1){
            if(any(levels(data[,rterms[select.terms]])!=levels(data[,rterms[1]]))){
              stop("terms involved in mm/str structures must have identical levels")
            }
            if(rtype=="mm"){
              colnames(Z[[k]])<-colnames(Z[[1]])
              Z[[1]]<-Z[[1]]+Z[[k]]
            }
            if(rtype=="str"){
              Z[[1]]<-cBind(Z[[1]],Z[[k]])
            }
          }
        }
        if(i==1){
          Zsave<-Z[[1]]
        }else{
          Zsave<-cBind(Zsave,Z[[1]])
        }
      }
      vnames<-colnames(X)  

      if(vtype=="us" | substr(vtype,1,3)=="cor" | substr(vtype,1,4)=="ante" | vtype=="sub"){
        if(Aterm==0 & rtype=="iid"){                     # if multivariate G-structure and no Ginv term drop levels where all effects are zero
          missing<-matrix(diff(Zsave@p)!=0, nrl, nfl)
          missing<-which(rowSums(missing)==0)
          if(length(missing)>0){
            Zsave<-Zsave[,-c(rep(missing, each=nfl)+rep((1:nfl-1)*nrl[1], length(missing))),drop=FALSE]
            nrl<-nrl-length(missing)
          }
        }
      }


      if(vtype=="idh"){
        Aterm<-rep(Aterm, nfl)
        vtype<-rep(vtype, nfl)
        nfl<-rep(1, nfl)
      }else{
        if(vtype[1]=="idv"){
          nrl<-nfl*nrl[1]
	  nfl<-1
	  vnames<-idv.vnames
        }
      }
      if(length(rterms)==0){
        rterms<-""
      }
      if(rtype=="iid"){
        if(vtype[1]=="us" | grepl("cor", vtype[1]) | substr(vtype[1],1,4)=="ante" | vtype[1]=="sub"){
          vnames<-paste(expand.grid(vnames, vnames)[,1],expand.grid(vnames, vnames)[,2], sep="MCMCsplit")
        }
        vnames<-paste(vnames, paste(rterms[select.terms], collapse=":"), sep=".")
      }
      if(rtype=="mm"){
        if(vtype[1]=="us" | grepl("cor", vtype[1]) | substr(vtype[1],1,4)=="ante" | vtype[1]=="sub"){
          vnames<-paste(expand.grid(vnames, vnames)[,1],expand.grid(vnames, vnames)[,2], sep="MCMCsplit")
        }
        vnames<-paste(vnames, paste(rterms, collapse="+"), sep=".")
      }
      if(rtype=="str"){
        nfl<-nfl*length(rterms)
        if(vtype[1]!="idh"){
          vnames<-paste(expand.grid(rterms, vnames)[,1],expand.grid(rterms, vnames)[,2], sep="MCMCsplit")
          vnames<-paste(expand.grid(vnames, vnames)[,1], expand.grid(vnames, vnames)[,2], sep=".")
        }else{
          vnamestmp<-vnames
          vnames<-c()
          for(i in 1:length(vnamestmp)){
            vnames<-c(vnames, paste(expand.grid(paste(rterms, vnamestmp[i], sep="MCMCsplit"),paste(rterms, vnamestmp[i], sep="MCMCsplit"))[,1], expand.grid(paste(rterms, vnamestmp[i], sep="MCMCsplit"),paste(rterms, vnamestmp[i], sep=":"))[,2], sep="."))
          }
        }
        vtype<-rep("us", length(vtype))
      }
      return(list(Z=Zsave, nfl=nfl, nrl=nrl, Aterm=Aterm, vtype=vtype, vnames=vnames, trait.ordering=trait.ordering))
    }
}
