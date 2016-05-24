makeZtX <-
  function(reinfo,Ku,Ju,uidx,Zmat=NULL){
    
    nKcol <- ncol(Ku)
    nJcol <- ncol(Ju)
    
    if(is.null(Zmat)){
      # one grouping variable
      
      numre <- length(reinfo$reNames[[1]])
      n2lev <- nlevels(reinfo$grpvar[,1])
      
      if(numre>1L){
        # multiple random effects
        
        ZtK <- matrix(0,n2lev*numre,nKcol)
        ZtJ <- matrix(0,n2lev*numre,nJcol)
        gsplt <- split(1:nrow(reinfo$grpvar),reinfo$grpvar[,1])
        recon <- 0L
        
        for(m in 1:numre){
          
          if(reinfo$reNames[[1]][m]=="1"){
            
            for(k in 1:n2lev){
              lidx <- length(gsplt[[k]])
              ZtK[k+recon,] <- colSums(matrix(Ku[uidx[gsplt[[k]]],],nrow=lidx,ncol=nKcol))
              ZtJ[k+recon,] <- colSums(matrix(Ju[uidx[gsplt[[k]]],],nrow=lidx,ncol=nJcol))
            }
            
          } else {
            
            cx <- match(reinfo$reNames[[1]][m],names(reinfo$raneff[[1]]))
            for(k in 1:n2lev){
              lidx <- length(gsplt[[k]])
              ZtK[k+recon,] <- colSums(reinfo$raneff[[1]][gsplt[[k]],cx]*matrix(Ku[uidx[gsplt[[k]]],],nrow=lidx,ncol=nKcol))
              ZtJ[k+recon,] <- colSums(reinfo$raneff[[1]][gsplt[[k]],cx]*matrix(Ju[uidx[gsplt[[k]]],],nrow=lidx,ncol=nJcol))
            }
            
          } # end if(reinfo$reNames[[1]][m]=="1")
          
          recon <- recon + n2lev
          
        } # end for(m in 1:numre)
        
        return(cbind(ZtK,ZtJ))
        
      } else {
        # single random effect
        
        ZtK <- matrix(0,n2lev,nKcol)
        ZtJ <- matrix(0,n2lev,nJcol)
        
        if(reinfo$reNames[[1]]=="1"){
          # random intercept only
          
          gsplt <- split(1:nrow(reinfo$grpvar),reinfo$grpvar[,1])
          for(k in 1:n2lev){
            lidx <- length(gsplt[[k]])
            ZtK[k,] <- colSums(matrix(Ku[uidx[gsplt[[k]]],],nrow=lidx,ncol=nKcol))
            ZtJ[k,] <- colSums(matrix(Ju[uidx[gsplt[[k]]],],nrow=lidx,ncol=nJcol))
          }
          
        } else if(is.factor(reinfo$raneff[[1]][,1])) {
          # subjects in groups
          
          gsplt <- split(1:nrow(reinfo$grpvar),list(reinfo$grpvar[,1],reinfo$raneff[[1]][,1]))
          glen <- sapply(gsplt, length)
          gsplt <- gsplt[glen>0L]
          for(k in 1:n2lev){
            lidx <- length(gsplt[[k]])
            ZtK[k,] <- colSums(matrix(Ku[uidx[gsplt[[k]]],],nrow=lidx,ncol=nKcol))
            ZtJ[k,] <- colSums(matrix(Ju[uidx[gsplt[[k]]],],nrow=lidx,ncol=nJcol))
          }
          
        } else {
          # random slope only
          
          gsplt <- split(1:nrow(reinfo$grpvar),reinfo$grpvar[,1])
          for(k in 1:n2lev){
            lidx <- length(gsplt[[k]])
            ZtK[k,] <- colSums(reinfo$raneff[[1]][gsplt[[k]],1]*matrix(Ku[uidx[gsplt[[k]]],],nrow=lidx,ncol=nKcol))
            ZtJ[k,] <- colSums(reinfo$raneff[[1]][gsplt[[k]],1]*matrix(Ju[uidx[gsplt[[k]]],],nrow=lidx,ncol=nJcol))
          }
          
        } # end if(reinfo$reNames[[1]]=="1")
        
        return(cbind(ZtK,ZtJ))
        
      } # end if(numre>1L)
      
    } else {
      return(crossprod(Zmat, cbind(Ku,Ju)[uidx,]))
    } # end if(is.null(Zmat))
    
    
  }