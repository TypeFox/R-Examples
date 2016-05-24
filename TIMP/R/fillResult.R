"fillResult" <- function (group, multimodel, thetalist, clpindepX, rlist,
	        rawtheta)
{
  irfvec <- rlist$irfvec
  QR.temp <- rlist$QR 
  psi <- rlist$psi 
  m <- multimodel@modellist 
  dset <- group[[1]][2]
  ## start getting the clp
  nnls <- multimodel@nnls 
  nnlscrit <- multimodel@nnlscrit 
  if(nnls) {
    X <- qr.X(QR.temp)
    if(length(nnlscrit$negpos) > 0) {
      con <- nnlscrit$spec[[as.character(group[[1]][1])]]
      cp <- coef(nnnpls(A = X, b = psi, con=con))
    }
    else {
      sol <- nnls(A = X, b = psi)
      cp <- coef(sol)
    }
    fitted <- X %*% cp 
    resid <- psi - fitted
  }  
  else {
    cp <- qr.coef(QR.temp, psi)
    resid <- qr.resid(QR.temp, psi)
    fitted <- qr.fitted(QR.temp, psi)   
  }
  cnt <- 1
  ## fill in results
  multimodel@nclp <- multimodel@nclp + length(cp)
  for (i in 1:length(group)) {
    dset <- group[[i]][2] 
    clpind <- group[[i]][1]
    res <- multimodel@fit@resultlist[[dset]]
    if(m[[dset]]@mod_type == "kin") { 
      irfvec[[i]][which(is.na(irfvec[[i]]))] <- 0
      if (length(m[[dset]]@cohspec$type) != 0) {
        if (m[[dset]]@cohspec$type ==  "freeirfdisp") 
          res@cohirf[[clpind]] <- rlist$cohirf[[i]] 
      }
      res@irfvec[[clpind]] <- irfvec[[i]] 
    }
    else
      res@irfvec[[clpind]] <- c(0,0) 
    if(m[[dset]]@clpType == "x2") 
      nt_or_nl <- m[[dset]]@nt
    else  nt_or_nl <- m[[dset]]@nl
    res@cp[[clpind]] <- cp
    res@fitted[[clpind]] <- fitted[cnt:(cnt + nt_or_nl - 1)]
    res@resid[[clpind]] <- resid[cnt:(cnt + nt_or_nl - 1)]        
    cnt <- cnt + nt_or_nl 
    multimodel@fit@rss <- multimodel@fit@rss + sum(res@resid[[clpind]]^2)
    multimodel@fit@resultlist[[dset]] <- res
  }
  multimodel
}

