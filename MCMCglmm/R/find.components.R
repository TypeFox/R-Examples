find.components<-function(x, data, nginverse=NULL){

  notrait=is.null(data$trait)

  vtype="us"    # form of covariances for the random effets of a random term classified by some factor 
  rtype="iid"   # form of covariances between random effets of random terms

  if(length(grep("^idh\\(", x))>0){
    vtype<-"idh"
  }
  if(length(grep("^cor\\(", x))>0){
    vtype<-"cor"
  }
  if(length(grep("^corh\\(", x))>0){
    vtype<-"corh"
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
  if(length(grep("(^|:)mm\\(", x))>0){
    rtype<-"mm"
  }
  if(length(grep("(^|:)str\\(", x))>0){
    rtype<-"str"
  }


  fformula<-gsub("^(us|cor|corh|corgh|corg|cors|idh|idv|ante.*)\\(", "", x)

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

  rterms<-strsplit(rterms, ":| \\+ ")[[1]]

  drop.rlevels<-rep(TRUE, length(rterms))

  fformula<-substr(fformula,1,closeB-1)
  fformula<-gsub(" ", "", fformula)

  interact<-strsplit(fformula, "\\+")[[1]]
  interact<-interact[grep(":|\\*",interact)]
  interact<-strsplit(interact, ":|\\*")

  if(fformula!=""){fformula<-c("+", fformula)}

  fformula<-as.formula(paste("~", "-1", paste(fformula, collapse=""), sep=""))

  ndata=data[1:min(dim(data)[1], 10),]

  if(notrait){
    ndata$trait<-gl(2, 1,min(dim(data)[1], 10))
  }

  X<-model.matrix(fformula, data=ndata)

  fformula=names(attr(X, "contrasts"))

  # IDEALLY - WHEN ANIMAL SPECIFIED at.level/set.level TERMS SHOULD BE RETAINED IN fformula

  if(is.null(fformula)==FALSE){
    if(any(fformula=="trait") & notrait){
      fformula<-fformula[-which(fformula=="trait")]
      if(length(fformula)==0){
        fformula=NULL
      }
    }
  }
  if(any(rterms%in%nginverse) | rtype!="iid"){
      drop.rlevels<-rep(FALSE, length(rterms))
  }

  if(is.null(fformula) | vtype=="idh"| vtype=="idv"){
    if(any(rterms%in%nginverse) | rtype!="iid"){
    }else{
      rterms<-NULL
    }
    fformula=NULL
  }

  # if any terms are interacted they have to be added
  if(length(interact)>0){
    i<-1
    while(i<=length(interact)){
      if(sum(interact[[i]]%in%fformula)>1){
        interact[[i]]<-interact[[i]][which(interact[[i]]%in%fformula)]
        fformula<-fformula[-which(fformula%in%interact[[i]])]
      }else{
        interact<-interact[-i]
        i<-i-1
      }
      i<-i+1
    }  
  }

  return(list(rterms=rterms,drop.rlevels=drop.rlevels, rtype=rtype,fformula=c(as.list(fformula), interact)))
}
