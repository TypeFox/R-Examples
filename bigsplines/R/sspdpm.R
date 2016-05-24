sspdpm <-
  function(xvars,type,rks,theknots,Etab,pred=FALSE,effect="all",intcpt=TRUE){
    ### make design and penalty matrix for bigssp
    
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
    for(j in 1:mfdim[2]){
      cidx <- which(Etab[,j]>0L)-1L
      lencidx <- length(cidx)
      if(lencidx==1L){
        if(type[cidx]!="prm"){
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
        } else if(effect=="all" || effect=="lin"){
          Kmat <- cbind(Kmat,xvars[[cidx]])
          Knames <- c(Knames,xnames[cidx])
        } # end if(type[chkidx]!="prm")
      } else if(lencidx==2L){
        pchk <- which(type[cidx]=="prm")
        lenpchk <- length(pchk)
        if(lenpchk==0L){
          if(effect=="all" || effect=="lin"){
            nkmt <- rkron(rks$kn[[cidx[1]]],rks$kn[[cidx[2]]])
            if(is.null(nkmt)==FALSE){
              Kmat <- cbind(Kmat,nkmt)
              Knames <- c(Knames,rep(paste(xnames[cidx[1]],xnames[cidx[2]],sep=":"),ncol(nkmt)))
            }
          }
          if(!is.null(rks$kn[[cidx[2]]])){
            if(effect=="all" || effect=="non"){
              Jmats <- cbind(Jmats,rks$jm[[cidx[1]]]*tcprod(rks$kn[[cidx[2]]],rks$qn[[cidx[2]]]))
              if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[1]]]*tcprod(rks$qn[[cidx[2]]]))}
              Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[2]],sep=":"))
            }
          }
          if(!is.null(rks$kn[[cidx[1]]])){
            if(effect=="all" || effect=="non"){
              Jmats <- cbind(Jmats,rks$jm[[cidx[2]]]*tcprod(rks$kn[[cidx[1]]],rks$qn[[cidx[1]]]))
              if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[2]]]*tcprod(rks$qn[[cidx[1]]]))}
              Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[2]],sep=":"))
            }
          }
          if(effect=="all" || effect=="non"){
            Jmats <- cbind(Jmats,rks$jm[[cidx[1]]]*rks$jm[[cidx[2]]])
            if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[1]]]*rks$qm[[cidx[2]]])}
            Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[2]],sep=":"))
          }
        } else if(lenpchk==1L){
          npchk <- (1:2)[-pchk]  
          if(!is.null(rks$kn[[cidx[npchk]]])){
            if(effect=="all" || effect=="lin"){
              Kmat <- cbind(Kmat,rks$kn[[cidx[npchk]]]*as.numeric(xvars[[cidx[pchk]]]))
              Knames <- c(Knames,rep(paste(xnames[cidx[1]],xnames[cidx[2]],sep=":"),xdim[cidx[npchk]]))
            }
          }
          if(effect=="all" || effect=="non"){
            Jmats <- cbind(Jmats,rks$jm[[cidx[npchk]]]*as.numeric(xvars[[cidx[pchk]]]))
            if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[npchk]]]*as.numeric(theknots[[cidx[pchk]]]))}
            Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[2]],sep=":"))
          }
        } else if(effect=="all" || effect=="lin"){
          Kmat <- cbind(Kmat,xvars[[cidx[1]]]*xvars[[cidx[2]]])
          Knames <- c(Knames,paste(xnames[cidx[1]],xnames[cidx[2]],sep=":"))
        } # end if(lenpchk==0L)
      } else {
        pchk <- which(type[cidx]=="prm")
        lenpchk <- length(pchk)
        if(lenpchk==0L){
          if(effect=="all" || effect=="lin"){
            nkmt <- rkron(rkron(rks$kn[[cidx[1]]],rks$kn[[cidx[2]]]),rks$kn[[cidx[3]]])
            if(!is.null(nkmt)){
              Kmat <- cbind(Kmat,nkmt)
              Knames <- c(Knames,rep(paste(xnames[cidx[1]],xnames[cidx[2]],xnames[cidx[3]],sep=":"),ncol(nkmt)))
            }
          }
          if(!is.null(rks$kn[[cidx[3]]])){
            if(!is.null(rks$kn[[cidx[2]]])){
              if(effect=="all" || effect=="non"){
                Jmats <- cbind(Jmats,rks$jm[[cidx[1]]]*tcprod(rks$kn[[cidx[2]]],rks$qn[[cidx[2]]])*tcprod(rks$kn[[cidx[3]]],rks$qn[[cidx[3]]]))
                if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[1]]]*tcprod(rks$qn[[cidx[2]]])*tcprod(rks$qn[[cidx[3]]]))}
                Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[2]],xnames[cidx[3]],sep=":"))
              }
            }
            if(!is.null(rks$kn[[cidx[1]]])){
              if(effect=="all" || effect=="non"){
                Jmats <- cbind(Jmats,rks$jm[[cidx[2]]]*tcprod(rks$kn[[cidx[1]]],rks$qn[[cidx[1]]])*tcprod(rks$kn[[cidx[3]]],rks$qn[[cidx[3]]]))
                if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[2]]]*tcprod(rks$qn[[cidx[1]]])*tcprod(rks$qn[[cidx[3]]]))}
                Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[2]],xnames[cidx[3]],sep=":"))
              }
            }
            if(effect=="all" || effect=="non"){
              Jmats <- cbind(Jmats,rks$jm[[cidx[1]]]*rks$jm[[cidx[2]]]*tcprod(rks$kn[[cidx[3]]],rks$qn[[cidx[3]]]))
              if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[1]]]*rks$qm[[cidx[2]]]*tcprod(rks$qn[[cidx[3]]]))}
              Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[2]],xnames[cidx[3]],sep=":"))
            }
          } # end if(is.null(rks$kn[[cidx[2]]])==FALSE)
          if(!is.null(rks$kn[[cidx[1]]])){
            if(!is.null(rks$kn[[cidx[2]]])){
              if(effect=="all" || effect=="non"){
                Jmats <- cbind(Jmats,rks$jm[[cidx[3]]]*tcprod(rks$kn[[cidx[1]]],rks$qn[[cidx[1]]])*tcprod(rks$kn[[cidx[2]]],rks$qn[[cidx[2]]]))
                if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[3]]]*tcprod(rks$qn[[cidx[1]]])*tcprod(rks$qn[[cidx[2]]]))}
                Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[2]],xnames[cidx[3]],sep=":"))
              }
            }
            if(effect=="all" || effect=="non"){
              Jmats <- cbind(Jmats,rks$jm[[cidx[2]]]*rks$jm[[cidx[3]]]*tcprod(rks$kn[[cidx[1]]],rks$qn[[cidx[1]]]))
              if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[2]]]*rks$qm[[cidx[3]]]*tcprod(rks$qn[[cidx[1]]]))}
              Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[2]],xnames[cidx[3]],sep=":"))
            }
          } # end if(is.null(rks$kn[[cidx[1]]])==FALSE)
          if(!is.null(rks$kn[[cidx[2]]])){
            if(effect=="all" || effect=="non"){
              Jmats <- cbind(Jmats,rks$jm[[cidx[1]]]*rks$jm[[cidx[3]]]*tcprod(rks$kn[[cidx[2]]],rks$qn[[cidx[2]]]))
              if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[1]]]*rks$qm[[cidx[3]]]*tcprod(rks$qn[[cidx[2]]]))}
              Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[2]],xnames[cidx[3]],sep=":"))
            }
          }
          if(effect=="all" || effect=="non"){
            Jmats <- cbind(Jmats,rks$jm[[cidx[1]]]*rks$jm[[cidx[2]]]*rks$jm[[cidx[3]]])
            if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[1]]]*rks$qm[[cidx[2]]]*rks$qm[[cidx[3]]])}
            Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[2]],xnames[cidx[3]],sep=":"))
          }
        } else if(lenpchk==1L){
          npchk <- (1:3)[-pchk]
          if(!any(c(is.null(rks$kn[[cidx[npchk[1]]]]),is.null(rks$kn[[cidx[npchk[2]]]])))){
            if(effect=="all" || effect=="lin"){
              nkmt <- rkron(rks$kn[[cidx[npchk[1]]]],rks$kn[[cidx[npchk[2]]]])*as.numeric(xvars[[cidx[pchk]]])
              Kmat <- cbind(Kmat,nkmt)
              Knames <- c(Knames,rep(paste(xnames[cidx[1]],xnames[cidx[2]],xnames[cidx[3]],sep=":"),ncol(nkmt)))
            }
          }
          if(is.null(rks$kn[[cidx[npchk[2]]]])==FALSE){
            if(effect=="all" || effect=="non"){
              Jmats <- cbind(Jmats,rks$jm[[cidx[npchk[1]]]]*tcprod(rks$kn[[cidx[npchk[2]]]],rks$qn[[cidx[npchk[2]]]])*as.numeric(xvars[[cidx[pchk]]]))
              if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[npchk[1]]]]*tcprod(rks$qn[[cidx[npchk[2]]]])*as.numeric(theknots[[cidx[pchk]]]))}
              Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[2]],xnames[cidx[3]],sep=":"))
            }
          }
          if(!is.null(rks$kn[[cidx[npchk[1]]]])){
            if(effect=="all" || effect=="non"){
              Jmats <- cbind(Jmats,rks$jm[[cidx[npchk[2]]]]*tcprod(rks$kn[[cidx[npchk[1]]]],rks$qn[[cidx[npchk[1]]]])*as.numeric(xvars[[cidx[pchk]]]))
              if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[npchk[2]]]]*tcprod(rks$qn[[cidx[npchk[1]]]])*as.numeric(theknots[[cidx[pchk]]]))}
              Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[2]],xnames[cidx[3]],sep=":"))
            }
          }
          if(effect=="all" || effect=="non"){
            Jmats <- cbind(Jmats,rks$jm[[cidx[npchk[1]]]]*rks$jm[[cidx[npchk[2]]]]*as.numeric(xvars[[cidx[pchk]]]))
            if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[npchk[1]]]]*rks$qm[[cidx[npchk[2]]]]*as.numeric(theknots[[cidx[pchk]]]))}
            Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[2]],xnames[cidx[3]],sep=":"))
          }
        } else if(lenpchk==2L){
          npchk <- (1:3)[-pchk]
          if(!is.null(rks$kn[[cidx[npchk]]])){
            if(effect=="all" || effect=="lin"){
              Kmat <- cbind(Kmat,rks$kn[[cidx[npchk]]]*as.numeric(xvars[[cidx[pchk[1]]]]*xvars[[cidx[pchk[2]]]]))
              Knames <- c(Knames,rep(paste(xnames[cidx[1]],xnames[cidx[2]],xnames[cidx[3]],sep=":"),xdim[cidx[npchk]]))
            }
          }
          if(effect=="all" || effect=="non"){
            Jmats <- cbind(Jmats,rks$jm[[cidx[npchk]]]*as.numeric(xvars[[cidx[pchk[1]]]]*xvars[[cidx[pchk[2]]]]))
            if(!pred){Qmats <- cbind(Qmats,rks$qm[[cidx[npchk]]]*as.numeric(theknots[[cidx[pchk[1]]]]*theknots[[cidx[pchk[2]]]]))}
            Jnames <- c(Jnames,paste(xnames[cidx[1]],xnames[cidx[2]],xnames[cidx[3]],sep=":"))
          }
        } else if(effect=="all" || effect=="lin"){
          Kmat <- cbind(Kmat,xvars[[cidx[1]]]*xvars[[cidx[2]]]*xvars[[cidx[3]]])
          Knames <- c(Knames,paste(xnames[cidx[1]],xnames[cidx[2]],xnames[cidx[3]],sep=":"))
        } # end if(lenpchk==0L)
      } # end if(lencidx==1L)
    } # end for(j in 1:nterms)
    
    if(!is.null(Kmat)){colnames(Kmat) <- Knames}
    if(!is.null(Jmats)){colnames(Jmats) <- rep(Jnames,each=nknots)}
    return(list(Kmat=Kmat,Jmats=Jmats,Qmats=Qmats))
    
  }