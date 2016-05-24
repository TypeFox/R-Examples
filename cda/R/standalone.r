.norm <- function(x) as.numeric(sqrt(t.default(x) %*% x))

.incident_field <- function(E0, k, r, angles){
  
  Ei <- matrix(0, ncol=nrow(angles), nrow=3*nrow(r))
  for (ii in seq(1, nrow(angles))){
    
    rot <- .euler(angles[ii,])
    k_r  <- t(rot) %*% k
    E0_r  <- t(rot) %*% E0
    kr <- r %*% k_r
    expikr <- exp(1i*kr)
    Ei[,ii] <- as.vector(matrix(c(E0_r[1]*expikr, E0_r[2]*expikr, 
                       E0_r[3]*expikr), nrow=3, byrow=T))
    
  }
  
  Ei
}

.extinction <- function(kn, P, Ei){
  N <- ncol(P)
  cext <- rep(0, N)
  
  for (ii in seq.int(N))
    cext[ii] <- Im(crossprod(Conj(Ei[,ii]), P[,ii]))
  
  4*pi*kn*cext
}

.absorption<- function(kn, P, Alpha){
  N <- ncol(P)
  cabs <- rep(0, N)
  Eexc <- Alpha %*% P
  
  for (ii in seq.int(N))
    cabs[ii] <- Im(crossprod(Conj(Eexc[,ii]), P[,ii])) - 
      kn^3*2/3* Re(crossprod(Conj(P[,ii]), P[,ii]))
  
  4*pi*kn*cabs
}

.euler <- function(x){
  alpha <- x[1]
  beta <- x[2]
  gamma <- x[3]
  Ra <- matrix(c(cos(alpha), sin(alpha), 0, -sin(alpha), cos(alpha), 0, 0, 0, 1), ncol=3, byrow=T)
  Rb <- matrix(c(1, 0, 0, 0, cos(beta), sin(beta), 0, -sin(beta), cos(beta)), ncol=3, byrow=T)
  Rc <- matrix(c(cos(gamma), sin(gamma), 0, -sin(gamma), cos(gamma), 0, 0, 0, 1), ncol=3, byrow=T)
  
  return(Rc%*%Rb%*%Ra)
  
}


.interaction_matrix <-  function(r, kn, beta, euler){
  
  N <- nrow(r)
  A <- matrix(0, 3*N, 3*N)
  
  I <- diag(3)
  
  for (jj in 1:N){
    for (kk in 1:N){
      
      if (jj != kk){          
        
        ind1 <- ((jj-1)*3+1):(jj*3)
        ind2 <- ((kk-1)*3+1):(kk*3)
        
        rk_to_rj <- r[jj,] - r[kk,]
        rjk <- .norm(rk_to_rj)
        rjk_hat <- rk_to_rj / rjk
        rjkrjk <- outer(rjk_hat, rjk_hat)
        
        Ajk <- exp(1i*kn*rjk) / rjk * (kn^2*(rjkrjk - I) + (1i*kn*rjk-1) / rjk^2 * (3*rjkrjk - I))           
        A[ind1, ind2] <- Ajk
        
      } else {
        
        ind1 <- ((jj-1)*3+1):(jj*3)
        rot <- .euler(euler[jj,])
        A[ind1, ind1]  <- t(rot) %*% diag(beta[ind1]) %*% rot
        
      }
    }
  }
  invisible(A)
}

.block_diagonal <-  function(beta, euler){
  
  N <- nrow(euler)
  A <- matrix(0, 3*N, 3*N)
  
  
  for (jj in 1:N){
    for (kk in 1:N){
      
      if (jj == kk){          
        
        ind1 <- ((jj-1)*3+1):(jj*3)
        rot <- .euler(euler[jj,])
        A[ind1, ind1]  <- t(rot) %*% diag(beta[ind1]) %*% rot
        
      }
    }
  }
  invisible(A)
}

