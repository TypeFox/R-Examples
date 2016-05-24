"getResidRet" <- 
  function(X, psi, rlist, returnX, finished, nnls, algorithm = "nls",
           nnlscrit=list(), group=0)
  {
    if(returnX)  return(as.vector(X))
    if(finished) {
      rlist$QR <- qr(X)
      rlist$psi <- psi 
      return(rlist) 
    }
    if(!nnls) { ## just varpro
      qty.temp <- qr.qty( qr(X) , psi)
      residQspace <- qty.temp[-(1:ncol(X))]
      retval <- residQspace
    }
    else {
      if(length(nnlscrit$negpos) > 0) {
        con <- nnlscrit$spec[[as.character(group[[1]][1])]]
        cp <- coef(nnnpls(A = X, b = psi, con=con))
      }
      else {
        sol <- try(nnls(A = X, b = psi))
        if(class(sol) == "try-error")
          cp <- rep(0, ncol(X))
        else
          cp <- coef(sol)
      }
      if(algorithm != "optim") 
        retval <- psi - X %*% cp
      else
        retval <-  sum((psi - X %*% cp)^2)
    }
    retval
  }
