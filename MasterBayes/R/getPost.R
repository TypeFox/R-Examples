"getPost"<-function(post, sP, X.list, nitt, thin, burnin, write_postG, write_postP, write_postA, unique_id, marker.type="MSW"){

########### extract posterior samples from C++ ######################################################################

  if(length(post$E1)!=0 & sP$estE1==TRUE){
    post$E1<-mcmc(t(matrix(post$E1, length(sP$E1), ceiling((nitt-burnin)/thin))))
    if(length(sP$E1)==1){
      colnames(post$E1)<-"E1"
    }else{
      colnames(post$E1)<-paste("E1.", names(sP$E1), sep="")
    }  
  }
  if(length(post$E2)!=0 & sP$estE2==TRUE){
    post$E2<-mcmc(t(matrix(post$E2, length(sP$E2), ceiling((nitt-burnin)/thin))))
    if(length(sP$E2)==1){
      colnames(post$E2)<-"E2"
    }else{
      colnames(post$E2)<-paste("E2.", names(sP$E2), sep="")
    }
  }
  
  if(length(post$beta)!=0 & sP$estbeta==TRUE){
    post$beta<-mcmc(t(matrix(post$beta, length(sP$beta), ceiling((nitt-burnin)/thin))))
    cnames<-unlist(lapply(X.list$X[[1]][seq(4,14,by=2)], colnames))
    Slinked<-c(grep("linked", cnames))
    if(length(Slinked)>0){
      cnames<-cnames[-Slinked[which(substr(names(cnames[Slinked]),2,2)=="S")]]
    }
    cnames<-paste(cnames, substr(names(cnames),2,3), sep=".")
    DSvar<-grep("S", substr(cnames, nchar(cnames), nchar(cnames)), ignore.case=FALSE)
    cnames<-substr(cnames,1, nchar(cnames)-1)
    cnames<-sapply(cnames, function(x){paste(strsplit(x, ".linked.D")[[1]], collapse="")})
    cnames<-sapply(cnames, function(x){paste(strsplit(x, "linked.")[[1]], collapse="")})
    if(length(DSvar)>0){
      cnames[DSvar]<-paste(cnames[DSvar], "S", sep="") 
    }
    colnames(post$beta)<-cnames
  }

  if(length(post$A)!=0 & sP$estA==TRUE & write_postA==TRUE){
    nall<-unlist(lapply(sP$A, length))
    count<-0
    post_AP<-as.list(1:length(sP$A))
    for(loc in 1:length(sP$A)){
      post_AP[[loc]]<-t(matrix(NA, nall[loc], ceiling((nitt-burnin)/thin)))
      colnames(post_AP[[loc]])<-names(sP$A[[loc]])        
      for(a_l in 1:nall[loc]){
        count<-count+1
        post_AP[[loc]][,a_l]<-post$A[seq(count,length(post$A),by=sum(nall))]
      }
    }
    post$A<-post_AP
    names(post$A)<-names(sP$A)
  }


  if(length(post$G)!=0 & sP$estG==TRUE & write_postG==TRUE){
    end<-0
    nall<-unlist(lapply(sP$A, length))
    nind<-length(unique_id)
    ngen<-(nall*(nall+1))/2
    post_GP<-as.list(1:length(sP$A))
    for(loc in 1:length(sP$A)){
      start<-end+1
      end<-end+(ngen[loc]*nind)
      post_GP[[loc]]<-matrix(post$G[start:end], nind, ngen[loc])
      if(marker.type=="MSC" || marker.type=="MSW"){
        gnames<-outer(names(sP$A[[loc]]), names(sP$A[[loc]]), paste, sep="/")
        colnames(post_GP[[loc]])<-c(t(gnames)[lower.tri(gnames, diag=TRUE)])
      }else{
        colnames(post_GP[[loc]])<-c("0/0", "0/1", "1/1")
      }
      rownames(post_GP[[loc]])<-unique_id
    }
    post$G<-post_GP
    names(post$G)<-names(sP$A)
  }


  if(length(post$P)!=0 & sP$estP==TRUE & write_postP=="JOINT"){
    post$P<-unique_id[post$P+1]
    post$P[is.na(post$P)==T]<-"us"
    post$P<-matrix(post$P, length(X.list$X), 2*ceiling((nitt-burnin)/thin))
    rownames(post$P)<-unique_id[as.numeric(names(X.list$X))]
  }

  if(length(post$P)!=0 & sP$estP==TRUE & write_postP=="MARGINAL"){
    Pmarginal<-as.list(1:length(X.list$X))
    dam_id<-unlist(lapply(X.list$X, function(x){x$restdam.id}))
    sire_id<-unlist(lapply(X.list$X, function(x){x$restsire.id}))
    ndam<-unlist(lapply(X.list$X, function(x){length(x$restdam.id)}))
    nsire<-unlist(lapply(X.list$X, function(x){length(x$restsire.id)}))
    start<-0
    starts<-0
    startd<-0
    for(i in 1:length(X.list$X)){ 
      Pmarginal[[i]]<-matrix(post$P[(1:(ndam[i]*nsire[i]))+start], nsire[i], ndam[i])
      colnames(Pmarginal[[i]])<-c(unique_id[dam_id[(1:ndam[i])+startd]])
      rownames(Pmarginal[[i]])<-c(unique_id[sire_id[(1:nsire[i])+starts]])
      startd<-startd+ndam[i]
      starts<-starts+nsire[i]
      start<-start+(ndam[i]*nsire[i])
      colnames(Pmarginal[[i]])[which(is.na(colnames(Pmarginal[[i]]))==T)]<-"us"
      rownames(Pmarginal[[i]])[which(is.na(rownames(Pmarginal[[i]]))==T)]<-"us"      
    }
    post$P<-Pmarginal
    names(post$P)<-unique_id[as.numeric(names(X.list$X))]
  }


########### create empty posterior samples for C++ ########################################################

  if(length(post$E1)==0 & sP$estE1==TRUE){
    post$E1<-rep(0, ceiling((nitt-burnin)/thin)*length(sP$E1))    # joint posterior distribution of E1 
  }
  if(length(post$E2)==0 & sP$estE2==TRUE){
    post$E2<-rep(0, ceiling((nitt-burnin)/thin)*length(sP$E2))    # joint posterior distribution of E1 
  }
  if(length(post$beta)==0 & sP$estbeta==TRUE){
    post$beta<-rep(0, ceiling((nitt-burnin)/thin)*length(sP$beta))
  }

  if(length(post$A)==0 & sP$estA==TRUE & write_postA==TRUE){
    post$A<-rep(0, ceiling((nitt-burnin)/thin)*length(unlist(sP$A)))
  }
  if(length(post$USdam)==0 & sP$estUSdam==TRUE){
    post$USdam<-rep(0, ceiling((nitt-burnin)/thin)*length(sP$USdam))
  }
  if(length(post$USsire)==0 & (sP$estUSsire==TRUE | sP$estUSsire=="USdam")){
    post$USsire<-rep(0, ceiling((nitt-burnin)/thin)*length(sP$USsire))
  }
  if(length(post$P)==0 & write_postP=="JOINT" & sP$estP==TRUE){
    post$P<-rep(0, ceiling((nitt-burnin)/thin)*length(X.list$X)*2)
  }
  if(length(post$P)==0 & write_postP=="MARGINAL" & sP$estP==TRUE){
    ndam<-unlist(lapply(X.list$X, function(x){length(x$restdam.id)}))
    nsire<-unlist(lapply(X.list$X, function(x){length(x$restsire.id)}))
    post$P<-rep(0, sum(ndam*nsire))
  }

  if(length(post$G)==0 & write_postG==TRUE & sP$estG==TRUE){
    nind<-length(unique_id)
    nall<-unlist(lapply(sP$A, length))
    ngen<-(nall*(nall+1))/2
    post$G<-rep(0, sum(nind*ngen)) 
  }
post
}
