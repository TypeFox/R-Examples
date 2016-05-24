ssadpm <-
  function(xvars,type,rks,theknots,Etab,pred=FALSE,effect="all",intcpt=TRUE){
    ### make design and penalty matrix for bigssa
    
    mfdim <- dim(Etab)
    nunewr <- NA
    nknots <- nrow(theknots[[1]])
    if(pred){nxvar <- length(xvars)} else{nxvar <- length(xvars)-1L}
    xnames <- rownames(Etab)[-1]
    xdim <- rep(NA,nxvar)
    for(k in 1:nxvar){
      if(!any(is.na(xvars[[k]]))){
        if(is.na(nunewr)){nunewr <- nrow(xvars[[k]])}
        xdim[k] <- ncol(xvars[[k]])
      }
    }
    if(intcpt){
      Kmat <- matrix(1,nunewr)
      Knames <- "0"
    } else{
      Kmat <- Knames <- NULL
    }
    Jmats <- Jnames <- Qmats <- NULL
    echk <- rowSums(Etab)[-1]
    
    for(j in 1:mfdim[2]){
      cidx <- which(Etab[,j]>0L)-1L
      lencidx <- length(cidx)
      if(lencidx==1L){
        if(!is.null(rks$qn[[cidx]])){
          if(effect=="all" || effect=="lin"){
            Kmat <- cbind(Kmat,rks$kn[[cidx]])
            Knames <- c(Knames,rep(xnames[cidx],xdim[cidx]))
          }
        }
        if(effect=="all" || effect=="non"){
          Jmats <- cbind(Jmats,rks$jm[[cidx]])
          Jnames <- c(Jnames,xnames[cidx])
        }
        if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx]])}
      } else if(lencidx==2L){
        if(effect=="all" || effect=="lin"){
          nkmt <- rkron(rks$kn[[cidx[1]]],rks$kn[[cidx[2]]])
          if(!is.null(nkmt)){
            Kmat <- cbind(Kmat,nkmt)
            Knames <- c(Knames,rep(paste(xnames[cidx[1]],xnames[cidx[2]],sep=":"),ncol(nkmt)))
          }
        }
        if(effect=="all" || effect=="non"){
          nidx <- match(xnames[cidx],Jnames)
          if(!is.null(rks$kn[[cidx[2]]])){
            if(is.na(nidx[1])){
              Jmats <- cbind(Jmats,rks$jm[[cidx[1]]]*tcprod(rks$kn[[cidx[2]]],rks$qn[[cidx[2]]]))
              Jnames <- c(Jnames,xnames[cidx[1]])
              if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[1]]]*tcprod(rks$qn[[cidx[2]]]))}
            } else {
              jidx <- (1+(nidx[1]-1)*nknots):(nidx[1]*nknots)
              Jmats[,jidx] <- Jmats[,jidx]+rks$jm[[cidx[1]]]*tcprod(rks$kn[[cidx[2]]],rks$qn[[cidx[2]]])
              if(!pred){Qmats[,jidx] <- Qmats[,jidx]+rks$qm[[cidx[1]]]*tcprod(rks$qn[[cidx[2]]])}
            }
          } # end if(is.null(rks$kn[[cidx[2]]])==FALSE)
          if(!is.null(rks$kn[[cidx[1]]])){
            if(is.na(nidx[2])){
              Jmats <- cbind(Jmats,rks$jm[[cidx[2]]]*tcprod(rks$kn[[cidx[1]]],rks$qn[[cidx[1]]]))
              Jnames <- c(Jnames,xnames[cidx[2]])
              if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[2]]]*tcprod(rks$qn[[cidx[1]]]))}
            } else {
              jidx <- (1+(nidx[2]-1)*nknots):(nidx[2]*nknots)
              Jmats[,jidx] <- Jmats[,jidx]+rks$jm[[cidx[2]]]*tcprod(rks$kn[[cidx[1]]],rks$qn[[cidx[1]]])
              if(!pred){Qmats[,jidx] <- Qmats[,jidx]+rks$qm[[cidx[2]]]*tcprod(rks$qn[[cidx[1]]])}
            }
          } # end if(is.null(rks$kn[[cidx[1]]])==FALSE)
          Jmats <- cbind(Jmats,rks$jm[[cidx[1]]]*rks$jm[[cidx[2]]])
          Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[2]],sep=":"))
          if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[1]]]*rks$qm[[cidx[2]]])}
        } # end if(effect=="all" || effect=="non")
        
      } else {
        
        if(effect=="all" || effect=="lin"){
          nkmt <- rkron(rkron(rks$kn[[cidx[1]]],rks$kn[[cidx[2]]]),rks$kn[[cidx[3]]])
          if(!is.null(nkmt)){
            Kmat <- cbind(Kmat,nkmt)
            Knames <- c(Knames,rep(paste(xnames[cidx[1]],xnames[cidx[2]],xnames[cidx[2]],sep=":"),ncol(nkmt)))
          }
        }
        
        if(effect=="all" || effect=="non"){
          nidx <- match(xnames[cidx],Jnames)
          inames <- c(paste(xnames[cidx[1]],xnames[cidx[2]],sep=":"),
                      paste(xnames[cidx[1]],xnames[cidx[3]],sep=":"),
                      paste(xnames[cidx[2]],xnames[cidx[3]],sep=":"))
          imatch <- match(inames,Jnames)
          if(!is.null(rks$kn[[cidx[3]]])){
            if(is.na(imatch[1])){
              Jmats <- cbind(Jmats,rks$jm[[cidx[1]]]*rks$jm[[cidx[2]]]*tcprod(rks$kn[[cidx[3]]],rks$qn[[cidx[3]]]))
              Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[2]],sep=":"))
              if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[1]]]*rks$qm[[cidx[2]]]*tcprod(rks$qn[[cidx[3]]]))}
            } else {
              jidx <- (1+(imatch[1]-1)*nknots):(imatch[1]*nknots)
              Jmats[,jidx] <- Jmats[,jidx]+rks$jm[[cidx[1]]]*rks$jm[[cidx[2]]]*tcprod(rks$kn[[cidx[3]]],rks$qn[[cidx[3]]])
              if(!pred){Qmats[,jidx] <- Qmats[,jidx]+rks$qm[[cidx[1]]]*rks$qm[[cidx[2]]]*tcprod(rks$qn[[cidx[3]]])}
            }
            if(!is.null(rks$kn[[cidx[2]]])){
              if(is.na(nidx[1])){
                Jmats <- cbind(Jmats,rks$jm[[cidx[1]]]*tcprod(rks$kn[[cidx[2]]],rks$qn[[cidx[2]]])*tcprod(rks$kn[[cidx[3]]],rks$qn[[cidx[3]]]))
                Jnames <- c(Jnames,xnames[cidx[1]])
                if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[1]]]*tcprod(rks$qn[[cidx[2]]])*tcprod(rks$qn[[cidx[3]]]))}
              } else {
                jidx <- (1+(nidx[1]-1)*nknots):(nidx[1]*nknots)
                Jmats[,jidx] <- Jmats[,jidx]+rks$jm[[cidx[1]]]*tcprod(rks$kn[[cidx[2]]],rks$qn[[cidx[2]]])*tcprod(rks$kn[[cidx[3]]],rks$qn[[cidx[3]]])
                if(!pred){Qmats[,jidx] <- Qmats[,jidx]+rks$qm[[cidx[1]]]*tcprod(rks$qn[[cidx[2]]])*tcprod(rks$qn[[cidx[3]]])}
              }
            } # end if(is.null(rks$kn[[cidx[2]]])==FALSE)
          } # end if(is.null(rks$kn[[cidx[3]]])==FALSE)
          
          if(!is.null(rks$kn[[cidx[2]]])){
            if(is.na(imatch[2])){
              Jmats <- cbind(Jmats,rks$jm[[cidx[1]]]*rks$jm[[cidx[3]]]*tcprod(rks$kn[[cidx[2]]],rks$qn[[cidx[2]]]))
              Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[3]],sep=":"))
              if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[1]]]*rks$qm[[cidx[3]]]*tcprod(rks$qn[[cidx[2]]]))}
            } else {
              jidx <- (1+(imatch[2]-1)*nknots):(imatch[2]*nknots)
              Jmats[,jidx] <- Jmats[,jidx]+rks$jm[[cidx[1]]]*rks$jm[[cidx[3]]]*tcprod(rks$kn[[cidx[2]]],rks$qn[[cidx[2]]])
              if(!pred){Qmats[,jidx] <- Qmats[,jidx]+rks$qm[[cidx[1]]]*rks$qm[[cidx[3]]]*tcprod(rks$qn[[cidx[2]]])}
            }
            if(!is.null(rks$kn[[cidx[1]]])){
              if(is.na(nidx[3])){
                Jmats <- cbind(Jmats,rks$jm[[cidx[3]]]*tcprod(rks$kn[[cidx[1]]],rks$qn[[cidx[1]]])*tcprod(rks$kn[[cidx[2]]],rks$qn[[cidx[2]]]))
                Jnames <- c(Jnames,xnames[cidx[3]])
                if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[3]]]*tcprod(rks$qn[[cidx[1]]])*tcprod(rks$qn[[cidx[2]]]))}
              } else {
                jidx <- (1+(nidx[3]-1)*nknots):(nidx[3]*nknots)
                Jmats[,jidx] <- Jmats[,jidx]+rks$jm[[cidx[3]]]*tcprod(rks$kn[[cidx[1]]],rks$qn[[cidx[1]]])*tcprod(rks$kn[[cidx[2]]],rks$qn[[cidx[2]]])
                if(!pred){Qmats[,jidx] <- Qmats[,jidx]+rks$qm[[cidx[3]]]*tcprod(rks$qn[[cidx[1]]])*tcprod(rks$qn[[cidx[2]]])}
              }
            } # end if(is.null(rks$kn[[cidx[1]]])==FALSE)
          } # end if(is.null(rks$kn[[cidx[2]]])==FALSE)
          
          if(!is.null(rks$kn[[cidx[1]]])){
            if(is.na(imatch[3])){
              Jmats <- cbind(Jmats,rks$jm[[cidx[2]]]*rks$jm[[cidx[3]]]*tcprod(rks$kn[[cidx[1]]],rks$qn[[cidx[1]]]))
              Jnames <- c(Jnames,paste(xnames[cidx[2]],xnames[cidx[3]],sep=":"))
              if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[2]]]*rks$qm[[cidx[3]]]*tcprod(rks$qn[[cidx[1]]]))}
            } else {
              jidx <- (1+(imatch[3]-1)*nknots):(imatch[3]*nknots)
              Jmats[,jidx] <- Jmats[,jidx]+rks$jm[[cidx[2]]]*rks$jm[[cidx[3]]]*tcprod(rks$kn[[cidx[1]]],rks$qn[[cidx[1]]])
              if(!pred){Qmats[,jidx] <- Qmats[,jidx]+rks$qm[[cidx[2]]]*rks$qm[[cidx[3]]]*tcprod(rks$qn[[cidx[1]]])}
            }
            if(!is.null(rks$kn[[cidx[3]]])){
              if(is.na(nidx[2])){
                Jmats <- cbind(Jmats,rks$jm[[cidx[2]]]*tcprod(rks$kn[[cidx[1]]],rks$qn[[cidx[1]]])*tcprod(rks$kn[[cidx[3]]],rks$qn[[cidx[3]]]))
                Jnames <- c(Jnames,xnames[cidx[2]])
                if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[2]]]*tcprod(rks$qn[[cidx[1]]])*tcprod(rks$qn[[cidx[3]]]))}
              } else {
                jidx <- (1+(nidx[2]-1)*nknots):(nidx[2]*nknots)
                Jmats[,jidx] <- Jmats[,jidx]+rks$jm[[cidx[2]]]*tcprod(rks$kn[[cidx[1]]],rks$qn[[cidx[1]]])*tcprod(rks$kn[[cidx[3]]],rks$qn[[cidx[3]]])
                if(!pred){Qmats[,jidx] <- Qmats[,jidx]+rks$qm[[cidx[2]]]*tcprod(rks$qn[[cidx[1]]])*tcprod(rks$qn[[cidx[3]]])}
              }
            } # end if(is.null(rks$kn[[cidx[1]]])==FALSE)
          } # end if(is.null(rks$kn[[cidx[2]]])==FALSE)
          
          Jmats <- cbind(Jmats,rks$jm[[cidx[1]]]*rks$jm[[cidx[2]]]*rks$jm[[cidx[3]]])
          Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[2]],xnames[cidx[3]],sep=":"))
          if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[1]]]*rks$qm[[cidx[2]]]*rks$qm[[cidx[3]]])}
          
        } # end if(effect=="all" || effect=="non")
        
      } # end if(lencidx==1L)
    } # end for(j in 1:mfdim[2])
    
    if(!is.null(Kmat)){colnames(Kmat) <- Knames}
    if(!is.null(Jmats)){colnames(Jmats) <- rep(Jnames,each=nknots)}
    return(list(Kmat=Kmat,Jmats=Jmats,Qmats=Qmats))
    
  }