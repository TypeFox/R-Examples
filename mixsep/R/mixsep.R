mixturesep <- function(x,nx=list(locus="sys",allele="allele",area="area",height="height"),m=2,
                       Jl=list(J2=J2,J20=J20,J3=J3,J30=J30),trace=TRUE,estimator=c("hat","tilde"),round=TRUE,
                       alternatives=TRUE,p=0.001,dropLocus=FALSE,fixedProfiles=list(),recur=FALSE,gui=FALSE,...){
  if(is.numeric(round)) round <- round
  else{ if(round) round <- 4 else round <- options()$digits }
  if(m>3) stop("The implementation can only handle up to three-person mixtures")
  if(m==2) J <- Jl[[1]]
  if(m==3) J <- Jl[[3]]
  y <- x
  est <- estimator[1]
  if(!is.element(est,c("hat","tilde"))) stop("Wrong type of estimator specified: 'estimator' needs to be either 'hat' or 'tilde'.")
  else{
    alfa <- get(paste("a",est,sep=""))
    tav <- get(paste("t",est,sep=""))
  }
  for(i in 1:length(nx)) names(x)[grep(nx[[i]],names(x))] <- names(nx)[i]
  x <- convertTab(x) ## converts all x columns to charaters
  x <- convertTab(x,mode="num",subset=c("area","height"))
  locusorder <- unique(paste(x$locus))
  x <- do.call("rbind",lapply(split(x,x$locus,drop=TRUE),function(y) within(y,{ns=nrow(y); asum=sum(area)})))
  if(length(ml <- unique(x$locus[x$ns>m*2]))>0)
    stop(paste("There are too many observed alleles (more than m*2=",m*2,") in locus ",paste(ml,sep=", ",collapse=", "),sep=""))
  if(length(fixedProfiles)>0){
    x <- handleFixed(x,f=fixedProfiles,m=m,gui=gui)
    if(is.null(x)) return(NULL)
    ## calling mixturesep itself in order to obtain best match with fewer known profiles
    if(recur){
      fm <- length(fixedProfiles)
      bm <- list()
      if(fm==1) bm <- c(bm,"U/U"=list(mixturesep(y,m=m,alternatives=FALSE,trace=FALSE,gui=gui)[c("profiles","stats")]))
      if(fm==2){
        bm <- c(bm,"F1/U"=list(mixturesep(y,m=m,alternatives=FALSE,trace=FALSE,gui=gui,fixedProfiles=list(fixedProfiles[[1]]))[c("profiles","stats")]))
        bm <- c(bm,"F2/U"=list(mixturesep(y,m=m,alternatives=FALSE,trace=FALSE,gui=gui,fixedProfiles=list(fixedProfiles[[2]]))[c("profiles","stats")]))
        bm <- c(bm,"U/U"=list(mixturesep(y,m=m,alternatives=FALSE,trace=FALSE,gui=gui)[c("profiles","stats")]))
      }
      if(m==3 & fm<3) names(bm) <- paste(names(bm),"U",sep="/")
      if(fm==3){
        bm <- c(bm,"F1/F2/U"=list(mixturesep(y,m=m,alternatives=FALSE,trace=FALSE,gui=gui,fixedProfiles=fixedProfiles[c(1,2)])[c("profiles","stats")]))
        bm <- c(bm,"F1/F3/U"=list(mixturesep(y,m=m,alternatives=FALSE,trace=FALSE,gui=gui,fixedProfiles=fixedProfiles[c(1,3)])[c("profiles","stats")]))
        bm <- c(bm,"F2/F3/U"=list(mixturesep(y,m=m,alternatives=FALSE,trace=FALSE,gui=gui,fixedProfiles=fixedProfiles[c(2,3)])[c("profiles","stats")]))
        bm <- c(bm,"F1/U/U"=list(mixturesep(y,m=m,alternatives=FALSE,trace=FALSE,gui=gui,fixedProfiles=list(fixedProfiles[[1]]))[c("profiles","stats")]))
        bm <- c(bm,"F2/U/U"=list(mixturesep(y,m=m,alternatives=FALSE,trace=FALSE,gui=gui,fixedProfiles=list(fixedProfiles[[2]]))[c("profiles","stats")]))
        bm <- c(bm,"F3/U/U"=list(mixturesep(y,m=m,alternatives=FALSE,trace=FALSE,gui=gui,fixedProfiles=list(fixedProfiles[[3]]))[c("profiles","stats")]))
        bm <- c(bm,"U/U/U"=list(mixturesep(y,m=m,alternatives=FALSE,trace=FALSE,gui=gui)[c("profiles","stats")]))
      }
    }
    Fixed <- grep("Fixed",names(x))
    LocusCol <- grep("locus",names(x))
    J <- Jl[[(m-1)*2]]
    x <- x[do.call("order",x[,c(LocusCol,Fixed)]),]
    x$Jmax <- NA
    x <- cbind(x,matrix(NA,nrow(x),m))
    names(x)[(ncol(x)-m+1):ncol(x)] <- paste("P",1:m,sep="")
    xs <- split(x,paste((2*m-x$ns),max(x$asum)-x$asum,x$locus,sep=":"))
    names(xs) <- unlist(lapply(strsplit(names(xs),":"),function(z) z[length(z)]))
  }
  else{
    x <- x[order(x$locus,x$area),]
    x$Jmax <- NA
    x <- cbind(x,matrix(NA,nrow(x),m))
    names(x)[(ncol(x)-m+1):ncol(x)] <- paste("P",1:m,sep="")
    xs <- split(x,paste((2*m-x$ns),max(x$asum)-x$asum,x$locus,sep=":"))
    names(xs) <- unlist(lapply(strsplit(names(xs),":"),function(z) z[length(z)]))
  }
  S <- length(xs)
  ## Iterations:
  L0 <- 0
  L1 <- -Inf
  iterations <- 0
  while(L0!=L1){
    if(trace) message(paste("Iteration:",format(iterations <- iterations+1,width=2),"\n"),appendLF=FALSE)
    L0 <- L1
    for(s in 1:S){
      locus.update <- locus.mixturesep(x=xs,s=s,J=J,m=m,est=est)
      xs[[s]] <- locus.update$xs
    }
    L1 <- locus.update$LIKE
  }
  if(dropLocus & length(fixedProfiles)==0){
    drop <- dropLoci(xs[match(locusorder,names(xs))],m=m,est=est,p=p)
    if(length(drop)==0) drop <- ""
    ## IMPLEMENT DROP-LOCUS FUNCTIONALITY
    ##     else{
    ##       xs <- xs[-(drop)]
    ##       locusorder <- locusorder[-drop]
    ##       for(s in 1:length(xs)) xs[[s]] <- locus.mixturesep(x=xs,s=s,J=J,m=m,est=est)$xs      
    ##     }
  }
  profiles <- lapply(xs,mmm)
  profiles <- profiles[match(locusorder,names(profiles))]
  profiles <- do.call("cbind",profiles)
  if(m==2) dimnames(profiles) <- list(c("Minor","Major"),locusorder)
  if(m==3) dimnames(profiles) <- list(c("Minor","Mid","Major"),locusorder)
  if(length(fixedProfiles)>0){
    dimnames(profiles)[[1]][1:length(grep("Fixed",names(x)))] <- names(x)[grep("Fixed",names(x))]
    if((mu <- m-length(fixedProfiles))>0)
      dimnames(profiles)[[1]][(length(grep("Fixed",names(x)))+1):nrow(profiles)] <- paste("Unknown",1:mu,sep="")
  }
  alpha <- alfa(xs,m=m)
  expAreas <- expectedAreas(xs,alpha=alpha,m=m)
  expAreas <- expAreas[order(ordered(expAreas$locus,levels=locusorder)),]
  if(length(fixedProfiles)==m) alternatives <- FALSE
  if(alternatives){
    alts <- alternative(xs,J=Jl[[(m-1)*2]],m=m,est=est,p=p)
    Fstats <- (lapply(alts,function(x) unlist(x$Fstats)))[match(locusorder,names(xs))]
    LocusTaus <- (lapply(alts,function(x) unlist(x$LocusTaus)))[match(locusorder,names(xs))]
    alts <- lapply(alts,function(x) x$alts)
    alts <- bindAlts(alts,m=m)
    noCombs <- alts$noCombs
    alts <- alts$alts[,match(locusorder,names(xs)),drop=FALSE]
    if(exists("bm")){
      BMstats <- as.data.frame(do.call("rbind",lapply(bm,function(z) z$stats)))
      BMstats$R2 <- round(min(BMstats$tau)/BMstats$tau,4)
      BMstats$R2[which.min(BMstats$tau)] <- ""
      BM <- do.call("rbind",lapply(bm,function(z) apply(z$profiles,2,paste,collapse="/")))
      if(all(alts=="")) alts <- BM
      else alts <- rbind(BM,alts)
    }
    dimnames(alts) <- list(alternatives=1:nrow(alts),Locus=locusorder)
    tauhat <- tav(xs,alpha=alpha,m=m)
    res <- list(profiles=profiles,alternatives=alts,noCombs=noCombs,Fstats=Fstats,LocusTaus=LocusTaus, ## was: data=y,
                expectedAreas=expAreas[,c("locus","allele","exp")],likelihood=L1,
                stats=round(c(alpha=alpha,tau=tauhat),round))
    if(exists("bm")){
      res$bm <- names(bm)
      res$bmstats <- BMstats
      res$R2 <- round(min(BMstats$tau)/tauhat,4)
    }
  }
  else{
    tauhat <- tav(xs,alpha=alpha,m=m)
    res <- list(profiles=profiles, ## was: data=y,
                expectedAreas=expAreas[,c("locus","allele","exp")],likelihood=L1,
                stats=round(c(alpha=alpha,tau=tauhat),round))
    res$alternatives <- matrix("",nrow=1,ncol=length(locusorder))
    res$noCombs <- NA
    if(exists("bm")){
      BMstats <- as.data.frame(do.call("rbind",lapply(bm,function(z) z$stats)))
      BMstats$R2 <- round(min(BMstats$tau)/BMstats$tau,4)
      BMstats$R2[which.min(BMstats$tau)] <- ""
      BM <- do.call("rbind",lapply(bm,function(z) apply(z$profiles,2,paste,collapse="/")))
#      BM <- do.call("rbind",lapply(bm,function(z) apply(z,2,paste,collapse="/")))
      dimnames(BM) <- list(BestMatch=1:nrow(BM),Locus=locusorder)
      res$bm <- names(bm)
      res$bmstats <- BMstats
      res$R2 <- round(min(BMstats$tau)/tauhat,4)
      res$alternatives <- BM
      res$noCombs <- 1
    }
   }
  if(dropLocus) res$dropLoci <- drop
  res
}

handleFixed <- function(x,f,m=2,gui=FALSE){
  fixed <- list()
  for(i in 1:length(f)){
    if(any(unlist(lapply(f[[i]],length))>2)){
      if(gui){
        tkmessageBox(message="A maximum of two alleles per profile is allowed for each locus",icon="error",type="ok")
        return(NULL)
      }
      else stop("A maximum of two alleles per profile is allowed for each locus")
    }
    tmpfixed <- data.frame(locus=rep(names(f[[i]]),unlist(lapply(f[[i]],length))),allele=unlist(f[[i]]),Fixed=1)
    tmpfixed <- tmpfixed[!duplicated(tmpfixed),]
    tmpfixed <- do.call("rbind",lapply(split(tmpfixed,tmpfixed$locus),function(y){ y$Fixed <- y$Fixed*2/sum(y$Fixed); y}))
    names(tmpfixed)[3] <- paste("Fixed",i,sep="")
    fixed <- c(fixed,list(tmpfixed))
  }
  for(i in 1:length(fixed)){
    x <- merge(x,fixed[[i]],by=c("locus","allele"),all.x=TRUE)
    x[is.na(x[,ncol(x)]),ncol(x)] <- 0
  }
  if(any(is.element(2,x[x$allele=="Y",grep("Fixed",names(x))]))){
    if(gui){
      tkmessageBox(message="It is not allowed to be homozygote YY in Amelogenin",icon="error",type="ok")
      return(NULL)
    }
    else stop("It is not allowed to be homozygote Y,Y in Amelogenin")
  }
  fixedSum <- unlist(lapply(split(x,x$locus),function(y) min(colSums(y[,grep("Fixed",names(y)),drop=FALSE]))))
  if(any(fixedSum==0)){
    if(gui){
      tkmessageBox(message=paste("Missing specification of alleles\nLocus: ",paste(names(fixedSum)[fixedSum==0],sep=", ",collapse=", "),sep=""),icon="error",type="ok")
      return(NULL)
    }
    else stop(paste("Missing specification of alleles\nLocus: ",paste(names(fixedSum)[fixedSum==0],sep=", ",collapse=", "),sep=""))
  }
  unknowns <- m-length(fixed)
  if(unknowns<0){
    if(gui){
      tkmessageBox(message="Too many fixed profiles specified compared to the number of contributors",icon="error",type="ok")
      return(NULL)
    }
    else stop("Too many fixed profiles specified compared to m")
  }
  else if(unknowns==0){
    fixedSum <- rowSums(x[,grep("Fixed",names(x)),drop=FALSE])
    if(any(fixedSum==0)){
      if(gui){
        tkmessageBox(message=paste("Not all alleles are accounted for by the fixed profiles\nLocus: ",paste(unique(x$locus[fixedSum==0]),collapse=", "),sep=""),icon="error",type="ok")
        return(NULL)
      }
      else stop(paste("Not all alleles are accounted for by the fixed profiles (locus: ",paste(unique(x$locus[fixedSum==0]),collapse=", "),")",sep=""))
    }
  }
  else{
    fixedSum <- unlist(lapply(split(x,x$locus),function(y) sum(rowSums(y[,grep("Fixed",names(y)),drop=FALSE])==0)))
    if(any(fixedSum>2*unknowns)){
      if(gui){
        tkmessageBox(message=paste("Too many alleles are unaccounted for by the fixed profiles with ",m," contributors\nLocus: ",paste(names(fixedSum)[fixedSum>2*unknowns],collapse=", "),sep=""),icon="error",type="ok")
        return(NULL)
      }
      else stop(paste("Too many alleles are unaccounted for by the fixed profiles with m=",m," (locus: ",paste(names(fixedSum)[fixedSum>2*unknowns],collapse=", "),")",sep=""))
    }
  }
  x
}

locus.mixturesep <- function(x,s,J,m,est){
  alfa <- get(paste("a",est,sep=""))
  tav <- get(paste("t",est,sep=""))
  PP <- paste("P",1:m,sep="")
  J <- J[[nrow(x[[s]])]]
  Fixed <- grep("Fixed",names(x[[s]]))
  if(length(Fixed)!=0 & nrow(x[[s]])>1){
    rmJ <- reduceJ(J,x[[s]])
    J <- J[-rmJ]
  }
  if(any(is.na(x[[s]]$P1))){
    x[[s]][,PP] <- J[[1]]
    x[[s]]$Jmax <- 1
    alpha <- alfa(x,m=m)
    TAU <- tav(x,alpha=alpha,m=m)
    if(m==3 & length(Fixed)==0){
      if( (alpha[1]>alpha[2]) | (alpha[1]>1-sum(alpha)) | (alpha[2]>1-sum(alpha)) ) TAU <- Inf
    }
    if(any(is.element(2,x[[s]][paste(x[[s]]$allele)=="Y",PP]))) TAU <- Inf ## Makes it impossible to be homozygote: Y,Y
    N <- sum(unlist(lapply(x,function(y) (!any(is.na(y$P1)))*(nrow(y)-1))))-1
    LIKE <- sqrt(TAU)^(-N) 
  }
  else{
    alpha <- alfa(x,m=m)
    TAU <- tav(x,alpha=alpha,m=m)
    if(m==3 & length(Fixed)==0){
      if( (alpha[1]>alpha[2]) | (alpha[1]>1-sum(alpha)) | (alpha[2]>1-sum(alpha)) ) TAU <- Inf
    }
    if(any(is.element(2,x[[s]][paste(x[[s]]$allele)=="Y",PP]))) TAU <- Inf ## Makes it impossible to be homozygote: Y,Y
    N <- sum(unlist(lapply(x,function(y) (!any(is.na(y$P1)))*(nrow(y)-1))))-1
    LIKE <- sqrt(TAU)^(-N)
  }
  jmax <- x[[s]]$Jmax[1]
  for(j in 1:length(J)){
    x[[s]][,PP] <- J[[j]]
    alpha <- alfa(x,m=m)
    tau <- tav(x,alpha=alpha,m=m)
    if(m==3 & length(Fixed)==0){
      if( (alpha[1]>alpha[2]) | (alpha[1]>1-sum(alpha)) | (alpha[2]>1-sum(alpha)) ) tau <- Inf
    }
    if(any(is.element(2,x[[s]][paste(x[[s]]$allele)=="Y",PP]))) tau <- Inf ## Makes it impossible to be homozygote: Y,Y
    like <- sqrt(tau)^(-N)
    if(like>LIKE){
      LIKE <- like
      jmax <- j
    }
  }
  x[[s]][,PP] <- J[[jmax]]
  x[[s]]$Jmax <- jmax
  list(xs=x[[s]],LIKE=LIKE)
}

reduceJ <- function(J,x){
  Fixed <- grep("Fixed",names(x))
  rmJ <- numeric()
  for(j in 1:length(J)){
    for(i in 1:length(Fixed)){
      if(any(J[[j]][,i]!=x[,Fixed[i]])){ rmJ <- c(rmJ,j); break }
    }
  }
  rmJ
}

bindAlts <- function(x,m=2){
  noCombs <- prod((unlist(lapply(x,length))+1))
  res <- lapply(x,function(a,m){
    if(length(a)==0) return(matrix("",1,1))
    if(length(a)==1) return(matrix(paste(mmm(a[[1]]),collapse="/"),1,1))
    else return(do.call("rbind",lapply(lapply(a,mmm),paste,collapse="/")))},m=m)
  res <- do.call("cbind",lapply(res,function(a,mnr){ rbind(a,matrix("",mnr-nrow(a),1))}, mnr=max(unlist(lapply(res,nrow)))))
  list(noCombs=noCombs,alts=res)
}

alternative <- function(x,J,m,est,p){
  alts <- lapply(x,function(...) NULL)
  for(s in 1:length(x)) alts[[s]] <- locus.alternative(x,s=s,J=J,m=m,est=est,p=p)
  alts
}

locus.alternative <- function(x,s,J,m,est,p=0.001){
  alts <- list()
  alts.sort <- list()
  Fstats <- numeric()
  Fstats.sort <- numeric()
  LocusTaus.sort <- numeric()
  alfa <- get(paste("a",est,sep=""))
  tav <- get(paste("t",est,sep=""))
  tav.locus <- get(paste("t",est,".locus",sep=""))
  PP <- paste("P",1:m,sep="")
  J <- J[[nrow(x[[s]])]]
  Fixed <- grep("Fixed",names(x[[s]]))
  if(length(Fixed)!=0){
    rmJ <- reduceJ(J,x[[s]])
    J <- J[-rmJ]
  }
  ## We already got a best match configuration in all loci
  jmax <- x[[s]]$Jmax[1]
  LocusTaus.min <- tav.locus(x[[s]],alpha=alfa(x,m=m),m=m)/(nrow(x[[s]])-1)
  ## Test for alternatives
  if(length(J)==1 | nrow(x[[s]])==1) return(list(alts=alts.sort,Fstats=Fstats,LocusTaus=c(LocusTaus.min,LocusTaus.sort)))
  else{
    for(j in (1:length(J))[-jmax]){
      x[[s]][,PP] <- J[[j]]
      alpha <- alfa(x,m=m)
      TAU <- tav(x,alpha=alpha,m=m)
      N <- sum(unlist(lapply(x,nrow)))-length(x)-1
      LIKE <- sqrt(TAU)^(-N)
      alphaNots <- alfa(x[-s],m=m)
      Ns <- N-(nrow(x[[s]])-1)
      TAUlocus <- sum(unlist(lapply(x[-s],tav.locus,alpha=alphaNots,m=m)))/Ns
      taulocus <- tav.locus(x[[s]],alpha=alphaNots,m=m)/(nrow(x[[s]])-1)
      if(m==3 & length(Fixed)==0){
        if( (alphaNots[1]>alphaNots[2]) | (alphaNots[1]>1-sum(alphaNots)) | (alphaNots[2]>1-sum(alphaNots)) ) taulocus <- Inf
      }
      if(any(is.element(2,x[[s]][paste(x[[s]]$allele)=="Y",PP]))) taulocus <- Inf ## Makes it impossible to be homozygote: Y,Y
      Fstat <- taulocus/TAUlocus
      Fstats <- c(Fstats,Fstat)
      alts <- c(alts,list(x[[s]]))
      if(Fstat < qf(p=1-p,df1=nrow(x[[s]])-1,df2=Ns)){
        Fstats.sort <- c(Fstats.sort,Fstat)
        LocusTaus.sort <- c(LocusTaus.sort, tav.locus(x[[s]],alpha=alfa(x,m=m),m=m)/(nrow(x[[s]])-1) )
        alts.sort <- c(alts.sort,list(x[[s]]))
      }
    }
    if(length(alts.sort)>0){
      alts.sort <- alts.sort[sort.list(Fstats.sort)]
      LocusTaus.sort <- LocusTaus.sort[sort.list(Fstats.sort)]
    }
  }
  list(alts=alts.sort,Fstats=Fstats,LocusTaus=c(LocusTaus.min,LocusTaus.sort))
}

dropLoci <- function(x,m,est,p=0.001){
  drop <- numeric()
  alfa <- get(paste("a",est,sep=""))
  tav <- get(paste("t",est,sep=""))
  tav.locus <- get(paste("t",est,".locus",sep=""))
  PP <- paste("P",1:m,sep="")
  S <- length(x)
  loclist <- 1:S
  while(length(loclist)>2){
    loclist <- loclist[!is.element(loclist,drop)]
    Fstats <- numeric()
    for(s in loclist){
      loclistNots <- loclist[loclist!=s]
      alphaNots <- alfa(x,m=m)
      N <- sum(unlist(lapply(x[loclistNots],nrow)))-length(loclist)-1
      Ns <- N-(nrow(x[[s]])-1)
      TAUlocus <- sum(unlist(lapply(x[loclistNots],tav.locus,alpha=alphaNots,m=m)))/Ns
      taulocus <- tav.locus(x[[s]],alpha=alphaNots,m=m)/(nrow(x[[s]])-1)
      Fstat <- taulocus/TAUlocus
      Fstats <- c(Fstats,Fstat)
    }
    cFs <- qf(p=1-p,df1=unlist(lapply(x[loclist],nrow))-1,df2=N-(unlist(lapply(x[loclist],nrow))-1))
    ## Bonferroni Correction of 'cFs' for multiple testing??
    if(any(Fstats>cFs)){
      drop <- c(drop,loclist[sort.list((Fstats-cFs), decreasing = TRUE)][1])
    }
    else return(drop)
  }
  (0:S)
}
