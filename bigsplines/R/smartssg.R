smartssg <-
  function(Etab,lambdas,family,Kmat,Jmats,Jnames,yvar,
           Qmats,nknots,ndpts,alpha,yty,nbf,fweights,fitweights,
           maxit,intol,subsamp,dispersion,gcvtype){
    
    ### get info
    rsums <- rowSums(Etab)[2:nrow(Etab)]
    if(any(rsums==0L)){
      ridx <- which(rsums==0L)
      Etab <- Etab[-(ridx+1),]
    }
    nxvar <- nrow(Etab)-1L
    gammas <- rep(0,nxvar)
    weights <- rep(0L,nxvar)
    xintidx <- rep(0L,nxvar)
    xnames <- rownames(Etab)[2:(nxvar+1)]
    thetas <- rep(NA,length(Jnames))
    
    ### check for 3-way interactions
    int3chk <- which(colSums(Etab)==3L)
    len3chk <- length(int3chk)
    if(len3chk>0){
      for(j in 1:len3chk){
        xchk <- which(Etab[,int3chk[j]]==1L)-1
        nchk <- xnames[xchk]
        n2way <- c(paste(nchk[1:2],collapse=":"),
                   paste(nchk[c(1,3)],collapse=":"),
                   paste(nchk[2:3],collapse=":"))
        n3way <- paste(nchk,collapse=":")
        jmidx <- c(match(n2way,Jnames),match(n3way,Jnames))
        for(k in 1:4){
          indx <- (1+(jmidx[k]-1)*nknots):(jmidx[k]*nknots)
          thetas[jmidx[k]] <- sum(diag(Qmats[,indx]))
        }
        gammas[xchk] <- gammas[xchk]+thetas[jmidx[3:1]]/thetas[jmidx[4]]
        weights[xchk] <- weights[xchk]+1L
        xintidx[xchk] <- xintidx[xchk]+1L
      } # end for(j in 1:len3chk)
    } # end if(len3chk>0)
    
    ### check for 2-way interactions
    int2chk <- which(colSums(Etab)==2L)
    len2chk <- length(int2chk)
    if(len2chk>0){
      for(j in 1:len2chk){
        xchk <- which(Etab[,int2chk[j]]==1L)-1
        if(any(xintidx[xchk]==0L)){
          nchk <- xnames[xchk]
          n2way <- paste(nchk,collapse=":")
          jmidx <- c(match(nchk,Jnames),match(n2way,Jnames))
          for(k in 1:3){
            if(is.na(thetas[jmidx[k]])){
              indx <- (1+(jmidx[k]-1)*nknots):(jmidx[k]*nknots)
              thetas[jmidx[k]] <- sum(diag(Qmats[,indx]))
            }
          }
          thenew <- thetas[jmidx[2:1]]/thetas[jmidx[3]]
          wchk <- which(xintidx[xchk]==0L)
          gammas[xchk[wchk]] <- gammas[xchk[wchk]]+thenew[wchk]
          weights[xchk[wchk]] <- weights[xchk[wchk]]+1L
        } # end if(any(xintidx[xchk]==0L))
      } # end for(j in 1:len2chk)
    } # end if(len2chk>0)
    
    ### check for empty 1-way
    if(any(weights==0L)){
      wchk <- which(weights==0L)
      lenw <- length(wchk)
      if(lenw==nxvar){
        # smart start (algorithm 3.2) from Gu and Wahba (1991)
        nfs <- dim(Qmats)[2]/nknots
        sdq <- rep(NA,nfs)
        for(jj in 1:nfs){
          jind <- ((jj-1)*nknots+1):(jj*nknots)
          sdq[jj] <- sum(diag(Qmats[,jind]))
        }
        gammas <- 1/sdq
        chat <- (lamcoefg(lambdas,gammas,family,Kmat,Jmats,yvar,
                          Qmats,nknots,ndpts,alpha,yty,nbf,fweights,
                          fitweights,maxit,intol,subsamp,dispersion,gcvtype))[[1]][(nbf+1):(nknots+nbf)]
        for(jj in 1:nfs){
          jind <- ((jj-1)*nknots+1):(jj*nknots)
          gammas[jj] <- (gammas[jj]^2)*crossprod(pdsXty(Qmats[,jind],chat))
        }
        return(gammas)
      } else {
        nzseq <- (1:nxvar)[-wchk]
        mgam <- mean(gammas[nzseq]/weights[nzseq])
        jmidx <- match(xnames[wchk],Jnames)
        for(k in 1:lenw){
          indx <- (1+(jmidx[k]-1)*nknots):(jmidx[k]*nknots)
          thetas[jmidx[k]] <- sum(diag(Qmats[,indx]))
        }
        gammas[wchk] <- mgam/thetas[jmidx]
        weights[wchk] <- 1
        
      } # end if(lenw==nxvar)
    } # end if(any(weights==0L))
    
    gammas <- gammas/weights
    
    ### fit fixed gamma model
    gamvec <- NULL
    for(j in 1:length(Jnames)){
      xi <- strsplit(Jnames[j],":")
      xidx <- match(xi[[1]],xnames)
      gamvec <- c(gamvec,prod(gammas[xidx]))
    }    
    chat <- (lamcoefg(lambdas,gamvec,family,Kmat,Jmats,yvar,
                      Qmats,nknots,ndpts,alpha,yty,nbf,fweights,
                      fitweights,maxit,intol,subsamp,dispersion,gcvtype))[[1]][(nbf+1):(nknots+nbf)]
    
    ### define vectors to collect new info
    gammasnew <- rep(0,nxvar)
    weightsnew <- rep(0L,nxvar)
    xintidxnew <- rep(0L,nxvar)
    thetasnew <- rep(NA,length(Jnames))
    
    ### check for 3-way interactions
    if(len3chk>0){
      for(j in 1:len3chk){
        xchk <- which(Etab[,int3chk[j]]==1L)-1
        nchk <- xnames[xchk]
        n2way <- c(paste(nchk[1:2],collapse=":"),
                   paste(nchk[c(1,3)],collapse=":"),
                   paste(nchk[2:3],collapse=":"))
        n3way <- paste(nchk,collapse=":")
        jmidx <- c(match(n2way,Jnames),match(n3way,Jnames))
        for(k in 1:4){
          indx <- (1+(jmidx[k]-1)*nknots):(jmidx[k]*nknots)
          thetasnew[jmidx[k]] <- (gamvec[jmidx[k]]^2)*crossprod(pdsXty(Qmats[,indx],chat))
        }
        gammasnew[xchk] <- gammasnew[xchk]+thetasnew[jmidx[4]]/thetasnew[jmidx[3:1]]
        weightsnew[xchk] <- weightsnew[xchk]+1L
        xintidxnew[xchk] <- xintidxnew[xchk]+1L
      } # end for(j in 1:len3chk)
    } # end if(len3chk>0)
    
    ### check for 2-way interactions
    if(len2chk>0){
      for(j in 1:len2chk){
        xchk <- which(Etab[,int2chk[j]]==1L)-1
        if(any(xintidxnew[xchk]==0L)){
          nchk <- xnames[xchk]
          n2way <- paste(nchk,collapse=":")
          jmidx <- c(match(nchk,Jnames),match(n2way,Jnames))
          for(k in 1:3){
            if(is.na(thetasnew[jmidx[k]])){
              indx <- (1+(jmidx[k]-1)*nknots):(jmidx[k]*nknots)
              thetasnew[jmidx[k]] <- (gamvec[jmidx[k]]^2)*crossprod(pdsXty(Qmats[,indx],chat))
            }
          }
          thenew <- thetasnew[jmidx[3]]/thetasnew[jmidx[2:1]]
          wchk <- which(xintidxnew[xchk]==0L)
          gammasnew[xchk[wchk]] <- gammasnew[xchk[wchk]]+thenew[wchk]
          weightsnew[xchk[wchk]] <- weightsnew[xchk[wchk]]+1L
        } # end if(any(xintidx[xchk]==0L))
      } # end for(j in 1:len2chk)
    } # end if(len2chk>0)
    
    ### check for empty 1-way
    if(any(weightsnew==0L)){
      wchk <- which(weightsnew==0L)
      nzseq <- (1:nxvar)[-wchk]
      jmidx <- match(xnames[wchk],Jnames)
      for(k in 1:length(wchk)){
        indx <- (1+(jmidx[k]-1)*nknots):(jmidx[k]*nknots)
        thetasnew[jmidx[k]] <- (gamvec[jmidx[k]]^2)*crossprod(pdsXty(Qmats[,indx],chat))
      }
      gammasnew[wchk] <- thetasnew[jmidx]
      weightsnew[wchk] <- 1
    } # end if(any(weights==0L))
    
    gammasnew/weightsnew
    
  }