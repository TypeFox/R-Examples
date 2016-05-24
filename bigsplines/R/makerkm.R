makerkm <-
  function(xvars,type,theknots,xrng,pred=FALSE,tpsinfo=NULL){
    ### make marginal reproducing kernel matrices
    
    nxvar <- length(xvars)
    if(pred){
      jm <- kn <- qn <- vector("list",nxvar)
      qm <- NA
    } else {
      jm <- kn <- qm <- qn <- vector("list",nxvar)
      tpsinfo <- vector("list",nxvar)
    }
    xdim <- rep(NA,nxvar)
    for(k in 1:nxvar){if(!any(is.na(xvars[[k]]))){xdim[k] <- ncol(xvars[[k]])}}
    nunewr <- NA
    nknots <- nrow(theknots[[1]])
    
    for(k in 1:nxvar){
      if(!is.na(xvars[[k]][[1]])){
        if(is.na(nunewr)){nunewr <- nrow(xvars[[k]])}
        if(type[[k]]=="cub"){
          kn[[k]] <- as.matrix(xvars[[k]]-0.5)
          jm[[k]] <- (.Fortran("cubker",xvars[[k]],theknots[[k]],nunewr,nknots,
                               matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[5]]
          qn[[k]] <- as.matrix(theknots[[k]]-0.5)
          if(pred==FALSE){
            qm[[k]] <- (.Fortran("cubkersym",theknots[[k]],nknots,
                                 matrix(0,nknots,nknots),PACKAGE="bigsplines"))[[3]]
          }
        } else if(type[[k]]=="per"){
          jm[[k]] <- (.Fortran("perker",xvars[[k]],theknots[[k]],nunewr,nknots,
                               matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[5]]
          if(pred==FALSE){
            qm[[k]] <- (.Fortran("perkersym",theknots[[k]],nknots,
                                 matrix(0,nknots,nknots),PACKAGE="bigsplines"))[[3]]
          }
        } else if(type[[k]]=="tps"){
          if(pred==FALSE){
            # get signing constant and projetion matrix
            if(xdim[k]<3){gconst <- 1} else{gconst <- (-1)}
            cqrd <- qr(cbind(1,theknots[[k]]),LAPACK=FALSE)
            CqrQ <- qr.Q(cqrd,complete=TRUE);     Rinv=solve(qr.R(cqrd))
            r1s <- sign(Rinv[1,1])
            if(r1s<0){
              CqrQ[,1] <- CqrQ[,1]*r1s
              Rinv[,1] <- Rinv[,1]*r1s
            }
            PhiQ <- sqrt(nknots)*CqrQ[,1:(xdim[k]+1)];   
            Pmat <- tcrossprod(CqrQ[,(xdim[k]+2):nknots])
            rm(CqrQ)
            # make penalty matrix
            Etilde <- gconst*(.Fortran("tpskersym",theknots[[k]],nknots,xdim[k],
                                       matrix(0,nknots,nknots),PACKAGE="bigsplines"))[[4]]
            qn[[k]] <- as.matrix(PhiQ[,2:(xdim[k]+1)])
            qm[[k]] <- Pmat%*%(Etilde%*%Pmat)
            # save tps info
            tpsinfo[[k]] <- list(Rinv,crossprod(PhiQ,Etilde),Pmat,qn[[k]])
            rm(Etilde)
            # make design matrix
            PhiX <- sqrt(nknots)*(cbind(1,xvars[[k]])%*%tpsinfo[[k]][[1]])
            kn[[k]] <- as.matrix(PhiX[,2:(xdim[k]+1)])
            jm[[k]] <- (.Fortran("tpsker",xvars[[k]],theknots[[k]],nunewr,
                                 xdim[k],nknots,matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[6]]
            jm[[k]] <- (gconst*jm[[k]]-((PhiX%*%tpsinfo[[k]][[2]])/nknots))%*%tpsinfo[[k]][[3]]
          } else {
            if(xdim[k]<3){gconst <- 1} else{gconst <- (-1)}
            PhiX <- sqrt(nknots)*(cbind(1,xvars[[k]])%*%tpsinfo[[k]][[1]])
            qn[[k]] <- tpsinfo[[k]][[4]]
            kn[[k]] <- as.matrix(PhiX[,2:(xdim[k]+1)])
            jm[[k]] <- (.Fortran("tpsker",xvars[[k]],theknots[[k]],nunewr,
                                 xdim[k],nknots,matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[6]]
            jm[[k]] <- (gconst*jm[[k]]-((PhiX%*%tpsinfo[[k]][[2]])/nknots))%*%tpsinfo[[k]][[3]]
          } # end if(pred==FALSE)
        } else if(type[[k]]=="ord"){
          jm[[k]] <- (.Fortran("ordker",xvars[[k]],theknots[[k]],nunewr,nknots,as.integer(xrng[[k]][2]),
                               matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[6]]
          if(pred==FALSE){
            qm[[k]] <- (.Fortran("ordkersym",theknots[[k]],nknots,as.integer(xrng[[k]][2]),
                                 matrix(0,nknots,nknots),PACKAGE="bigsplines"))[[4]]
          }
        } else if(type[[k]]=="nom"){
          jm[[k]] <- (.Fortran("nomker",xvars[[k]],theknots[[k]],nunewr,nknots,1/xrng[[k]][2],
                               matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[6]]
          if(pred==FALSE){
            qm[[k]] <- (.Fortran("nomkersym",theknots[[k]],nknots,1/xrng[[k]][2],
                                 matrix(0,nknots,nknots),PACKAGE="bigsplines"))[[4]]
          }
        } else if(type[[k]]=="cub0"){
          kn[[k]] <- as.matrix(xvars[[k]])
          jm[[k]] <- (.Fortran("cubkerz",xvars[[k]],theknots[[k]],nunewr,nknots,
                               matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[5]]
          qn[[k]] <- as.matrix(theknots[[k]])
          if(pred==FALSE){
            qm[[k]] <- (.Fortran("cubkerzsym",theknots[[k]],nknots,
                                 matrix(0,nknots,nknots),PACKAGE="bigsplines"))[[3]]
          }
        } # end if(type[[k]]=="cub")
      } # end if(is.na(xvars[[k]][[1]])==FALSE)
    } # end for(k in 1:nxvar)
    
    rks <- list(kn=kn,jm=jm,qn=qn,qm=qm,tpsinfo=tpsinfo)
    return(rks)
    
  }