mixedVarSim <-
  function (X1,X2) {  #cos2 of canonical analysis
    
    n <- length(X1)
    
    if ((is.numeric(X1) && (is.numeric(X2)))) { #cas quanti-quanti
      Z1 <- recodquant(X1)$Z
      Z2 <- recodquant(X2)$Z
      sim <- (t(Z1)%*%Z2/n)^2
    }
    if ((is.numeric(X1) && (!is.numeric(X2)))) { #cas quanti-quali
      Z1 <- recodquant(X1)$Z
      G2 <- recodqual(X2)
      ns <- apply(G2,2,sum)
      A <- t(G2)%*%Z1/ns
      sim <- sum((A^2*ns/n)) 
    }
    if ((!is.numeric(X1) && (is.numeric(X2)))) { #cas quali-quanli
      G1 <- recodqual(X1)
      Z2 <- recodquant(X2)$Z
      ns <- apply(G1,2,sum)
      A <- t(G1)%*%Z2/ns
      sim <- sum((A^2*ns/n)) 
    }
    if ((!is.numeric(X1) && (!is.numeric(X2)))) { #cas quali-quali
      G1 <- recodqual(X1)
      ns <- apply(G1,2,sum)
      ps <- ns/n
      X1 <- sweep(G1,MARGIN=2,STATS=sqrt(ps),FUN="/")
      G2 <- recodqual(X2)
      ns <- apply(G2,2,sum)
      ps <- ns/n
      X2 <- sweep(G2,MARGIN=2,STATS=sqrt(ps),FUN="/")
      r <- ncol(X1)
      s <- ncol(X2)
      m <- which.min(c(n,r,s))
      if (m==1) {
        A1 <- X1%*%t(X1)/n
        A2 <- X2%*%t(X2)/n
        A <- A1%*%A2
        e <- eigen(A)
        sim <- Re(e$values[2]) 
      } else {
        V12 <- t(X1)%*%X2/n
        V21 <- t(X2)%*%X1/n
        if (m==2) V<-V12%*%V21
        if (m==3) V<-V21%*%V12
        e <- eigen(V)
        sim <- Re(e$values[2]) 
      }	
    }		
    return(sim)
  }

