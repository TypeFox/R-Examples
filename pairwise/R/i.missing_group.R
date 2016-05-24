### ------------- missing_group ----------------------

missing_group<-function(X,s=NULL,all=FALSE){
  # aufteilung der daten X nach missing pattern X ist ein data.frame
  # sortierung der missing pattern gruppen (optional) nach summe der itemschwierigkeiten
  
  if(length(s)==0){s=rep(0,dim(X)[2])}
  s<-(s-min(s))+1
  e0<-apply(is.na(X)*1, 1, function(x){ paste(x, collapse="") })
  e1<-unique(e0)
  e2<-lapply(strsplit(e1,"",fixed = TRUE),as.numeric)
  ######################################
  reihenfolge<-order(unlist(lapply(e2,function(x){sum(x*s)})))
  ######################################
  mis_pattern<-unlist(lapply(lapply(e2,function(x){ifelse(x,"f","x")}), paste, collapse =""))[reihenfolge]
  person_index<-(lapply(as.list(e1),function(x){which(e0==x)}))[reihenfolge]
  names(person_index)<-mis_pattern
  mis_index=e2[reihenfolge]
  names(mis_index)<-mis_pattern
  use_index<-lapply(e2,function(x){which(x==0)})[reihenfolge]
  names(use_index)<-mis_pattern
  
  n<-unlist(lapply(person_index,length))
  size<-n/dim(X)[1]
  # daten X (alle items) als list nach missing pattern 
  X_mis_group<-mapply(function(x){as.matrix(X[x,])},person_index,SIMPLIFY = FALSE)
  
  X_nonmis_group<-mapply(function(x,y){as.matrix(X[x,y])},person_index,use_index,SIMPLIFY = FALSE)
  
  erg<-list(X.mis.group=X_mis_group,use.index=use_index,X.nonmis.group=X_nonmis_group,missing.index=mis_index, person.index=person_index,missing.pattern=mis_pattern,n=n,size.n=size)
  if(all==FALSE){ # rausenehmen der NA gruppen mit weniger als 2 items.
    ben <- names(which(lapply(erg$use.index,length)>=2)) #benutzen 
    nben <- which(lapply(erg$use.index,length)<2) #nicht benutzen
    erliste <- erg[1:5] # für listen elemente
    er1_ben <- list()
    er1_nben <- list()
    for (i in 1:5){er1_ben[[i]] <- erliste[[i]][ben] }
    for (i in 1:5){er1_nben[[i]] <- erliste[[i]][nben] }
    names(er1_ben) <- names(erg[1:5])
    names(er1_nben) <- names(erg[1:5])
    ervect <- erg[7:8]
    er2_ben <- lapply(ervect,function(x){x[ben]    })
    er2_nben <- lapply(ervect,function(x){x[nben]    })
    
    Pin <- er1_ben
    Pin[[6]] <- (ben)
    Pin[7:8] <- er2_ben
    names(Pin) <- names(erg)
    Pout <- er1_nben
    Pout[[6]] <- (nben)
    Pout[7:8] <- er2_nben
    names(Pout) <- names(erg)
    
    result<-list(Pin=Pin,Pout=Pout)
    erg<-result
  }
  return(erg)
}