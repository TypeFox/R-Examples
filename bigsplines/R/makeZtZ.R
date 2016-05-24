makeZtZ <-
  function(grpvar,reNames,raneff,yvar){
    
    ngvar <- ncol(grpvar)
    Zmat <- rdm <- znames <- NULL
    ndpts <- length(yvar)
    
    if(ngvar>1L){
      # multiple grouping variables
      
      # initialize Zmat
      ncolZ <- sum(sapply(reNames,length)*apply(grpvar,2,function(x) nlevels(factor(x))))
      Zmat <- matrix(0,ndpts,ncolZ)
      zidx <- 0
      
      # fill-in Zmat
      for(j in 1:ngvar){
        
        numre <- length(reNames[[j]])
        n2lev <- nlevels(grpvar[,j])
        gsplt <- split(1:ndpts,grpvar[,j])
        
        if(numre>1L){
          # multiple random effects
          
          for(m in 1:numre){
            
            zidx <- zidx[length(zidx)] + 1:n2lev
            if(reNames[[j]][m]=="1"){
              for(k in 1:n2lev) Zmat[gsplt[[k]],zidx[k]] <- 1
            } else {
              cx <- match(reNames[[j]][m],names(raneff[[j]]))
              for(k in 1:n2lev) Zmat[gsplt[[k]],zidx[k]] <- raneff[[j]][gsplt[[k]],cx]
            }
            znames <- c(znames, paste(levels(grpvar[,j]),reNames[[j]][m],sep="."))
            
          } # end for(m in 1:numre)
          
          newrdm <- rep(n2lev,numre)
          names(newrdm) <- paste(names(grpvar)[j],reNames[[j]],sep=".")
          rdm <- c(rdm,newrdm)
          
        } else {
          # single random effect
          
          zidx <- zidx[length(zidx)] + 1:n2lev
          
          if(reNames[[j]]=="1"){
            # random intercept only
            
            for(k in 1:n2lev) Zmat[gsplt[[k]],zidx[k]] <- 1
            znames <- c(znames, levels(grpvar[,j]))
            newrdm <- n2lev
            names(newrdm) <- paste(names(grpvar)[j],reNames[[j]],sep=".")
            rdm <- c(rdm,newrdm)
            
          } else if(is.factor(raneff[[j]][,1])) {
            # subjects in groups
            
            n1lev <- nlevels(raneff[[j]][,1])
            ysp <- split(yvar,list(raneff[[j]][,1],grpvar[,j]))
            zt <- matrix(unlist(lapply(ysp,length)),n1lev,n2lev)
            zt[zt==0] <- NA
            nax <- apply(zt,2,function(x) sum(is.na(x)))
            if(any(nax!=(nrow(zt)-1L))) stop(paste("Group 2 variable",names(grpvar)[j],"must be nested within Group 1 variable",reNames[[j]],"."))
            gix <- NULL
            for(k in 1:n1lev) gix <- c(gix, which(!is.na(zt[k,])))
            gsplt <- gsplt[gix]
            for(k in 1:n2lev) Zmat[gsplt[[k]],zidx[k]] <- 1
            znames <- c(znames, levels(grpvar[,j])[gix])
            newrdm <- rowSums(!is.na(zt))
            names(newrdm) <- paste(names(grpvar)[j],reNames[[j]],levels(raneff[[j]][,1]),sep=".")
            rdm <- c(rdm,newrdm)
            
          } else {
            # random slope only
            
            for(k in 1:n2lev) Zmat[gsplt[[k]],zidx[k]] <- raneff[[j]][gsplt[[k]],1]
            znames <- c(znames, levels(grpvar[,j]))
            newrdm <- n2lev
            names(newrdm) <- paste(names(grpvar)[j],reNames[[j]],sep=".")
            rdm <- c(rdm,newrdm)
            
          } # end if(reNames[[1]]=="1")
          
        } # end if(numre>1L)
        
      } # end for(j in 1:ngvar)
      
      ZtZ <- crossprod(Zmat)
      Zty <- crossprod(Zmat,yvar)
      names(Zty) <- znames
      
    } else {
      # single grouping variable
      
      numre <- length(reNames[[1]])
      n2lev <- nlevels(grpvar[,1])
      
      if(numre>1L){
        # multiple random effects
        
        # get ZtZ and Zty
        if(any(reNames[[1]]=="1")){
          ix <- which(reNames[[1]]=="1")
          if(ix!=1L) stop("Random intercept should precede random slopes when specifying random formula.\n Example:  use ~1+x|group instead of ~x+1|group.")
          ysp <- split(cbind(yvar,1,raneff[[1]]),grpvar[,1])
        } else {
          ysp <- split(cbind(yvar,raneff[[1]]),grpvar[,1])
        }
        ZtZlist <- lapply(ysp,function(x) crossprod(as.matrix(x[,-1])))
        ZtZ <- matrix(0,numre*n2lev,numre*n2lev)
        idx <- 1:numre
        for(m in 1:n2lev){
          ZtZ[idx,idx] <- ZtZlist[[m]]
          idx <- idx + numre
        }
        Zty <- unlist(lapply(ysp,function(x) crossprod(as.matrix(x[,-1]),x[,1])))
        
        # reorder for remlvc function
        zseq <- seq(1,n2lev*numre,by=numre)
        znames <- zidx <- NULL
        for(k in 1:numre) {
          znames <- c(znames, paste(levels(grpvar[,1]),reNames[[1]][k],sep="."))
          zidx <- c(zidx, zseq+k-1)
        }
        ZtZ <- ZtZ[zidx,zidx]
        Zty <- Zty[zidx]
        names(Zty) <- znames
        rdm <- rep(n2lev,numre)
        names(rdm) <- paste(names(grpvar),reNames[[1]],sep=".")
        
      } else {
        # single random effect
        
        if(reNames[[1]]=="1"){
          
          ysp <- split(yvar,grpvar[,1])
          ZtZ <- unlist(lapply(ysp,length))
          Zty <- unlist(lapply(ysp,sum))
          names(Zty) <- levels(grpvar[,1])
          rdm <- length(Zty)
          names(rdm) <- paste(names(grpvar),reNames[[1]],sep=".")
          
        } else if(is.factor(raneff[[1]][,1])) {
          # subjects in groups
          
          #browser(1>0)
          
          n1lev <- nlevels(raneff[[1]][,1])
          ysp <- split(yvar,list(raneff[[1]][,1],grpvar[,1]))
          zt <- matrix(unlist(lapply(ysp,length)),n1lev,n2lev)
          zt[zt==0] <- NA
          nax <- apply(zt,2,function(x) sum(is.na(x)))
          if(any(nax!=(nrow(zt)-1L))) stop(paste("Group 2 variable",names(grpvar),"must be nested within Group 1 variable",reNames[[1]],"."))
          gix <- NULL
          for(k in 1:n1lev) gix <- c(gix, which(!is.na(zt[k,])))
          ZtZ <- colSums(zt[,gix], na.rm=TRUE)
          Zty <- colSums(matrix(unlist(lapply(ysp,sum)),n1lev,n2lev))[gix]
          names(Zty) <- levels(grpvar[,1])[gix]
          rdm <- rowSums(!is.na(zt))
          names(rdm) <- paste(names(grpvar),reNames[[1]],levels(raneff[[1]][,1]),sep=".")
          
        } else {
          # random slope only
          
          ysp <- split(cbind(yvar,raneff[[1]]),grpvar[,1])
          ZtZ <- unlist(lapply(ysp,function(x) sum(x[,2]^2)))
          Zty <- unlist(lapply(ysp,function(x) sum(x[,1]*x[,2])))
          names(Zty) <- levels(grpvar[,1])
          rdm <- length(Zty)
          names(rdm) <- paste(names(grpvar),reNames[[1]],sep=".")
          
        } # end if(reNames[[1]]=="1")
        
      } # end if(length(reNames[[1]])>1L)
      
    } # end if(ncol(grpvar)>1L)
    
    return(list(Zty=Zty,ZtZ=ZtZ,rdm=rdm,Zmat=Zmat))
    
  }