profile.newton.se <-
function(x0, D1, D2, r, n1, n2, control = NULL) {
  if(is.null(control)){   
    lklMethod = 'REML'      
    lb = 0.01   	    	
    ub = 100 		       	
    s2profile = 'se' 		
    tol = 0.01          
  } else{    
    lklMethod = control$lklMethod
    lb = control$lb
    ub = control$ub
    s2profile = control$s2profile
    tol = control$tol
  }
  
  p = dim(D1)[1]
  Ip = diag(rep(1, p), p, p)
  n = n1 + n2
  
  ##This is based on the notation in Lindstrom and Bates (1988)
  N = n * p 
  
  ##residual matrices
  R1 = matrix(0, p, p)
  R2 = matrix(0, p, p)
  for (i in 1:n1) {
    R1 = R1 + r[, i] %o% r[, i]
  }
  for (i in (n1 + 1):n) {
    R2 = R2 + r[, i] %o% r[, i]
  }
  
  D1D1 = (D1 %*% t(D1))
  D2D2 = (D2 %*% t(D2))
  
  matTr <- function(z) sum(diag(z))
  
  #obj fn
  f <- function(tau) {
    V1 = tau * D1D1 + Ip
    V2 = tau * D2D2 + Ip
    V1inv = chol2inv(chol(V1))
    V2inv = chol2inv(chol(V2))
    
    tmp = ifelse((lklMethod == "REML"), N - 2*p, N)
    
    val = n1 * as.numeric(determinant(V1)$modulus) + n2 * as.numeric(determinant(V2)$modulus) + 
      tmp * log(matTr(V1inv %*% R1) + matTr(V2inv %*% R2))
    
    ##for REML
    if (lklMethod == "REML") {
      tmp = n1 * (t(D1) %*% V1inv %*% D1) 
      val = val + as.numeric(determinant(tmp)$modulus)
      tmp = n2 * (t(D2) %*% V2inv %*% D2)
      val = val + as.numeric(determinant(tmp)$modulus)
    }
    
    return(val)
  }
  
  #gradient fn
  g <- function(tau) {
    V1 = tau * D1D1 + Ip
    V2 = tau * D2D2 + Ip
    dV1 = D1D1
    dV2 = D2D2
    V1inv = chol2inv(chol(V1))
    V2inv = chol2inv(chol(V2))
    C1 = t(D1) %*% V1inv %*% D1
    C2 = t(D2) %*% V2inv %*% D2
    
    tmp = ifelse((lklMethod == "REML"), N - 2*p, N)
    
    val = n1 * matTr(C1) + n2 * matTr(C2) - 
      tmp * (matTr(V1inv %*% dV1 %*% V1inv %*% R1) + matTr(V2inv %*% dV2 %*% V2inv %*% R2))/(matTr(V1inv %*% R1) + matTr(V2inv %*% R2))
    ##for REML
    if (lklMethod == "REML") {
      val = val - matTr(C1) - matTr(C2) 
    }
    
    return(val)
  }
  
  #hessian fn
  
  h <- function(tau) {
    V1 = tau * D1D1 + Ip
    V2 = tau * D2D2 + Ip
    dV1 = D1D1
    dV2 = D2D2
    V1inv = chol2inv(chol(V1))
    V2inv = chol2inv(chol(V2))
    C1 = t(D1) %*% V1inv %*% D1
    C2 = t(D2) %*% V2inv %*% D2
    C1C1 = C1 %*% C1
    C2C2 = C2 %*% C2
    
    tmp = ifelse((lklMethod == "REML"), N - 2*p, N)
    
    hes1 = (matTr(V1inv %*% D1 %*% C1 %*% t(D1) %*% V1inv %*% R1) + matTr(V2inv %*% D2 %*% C2 %*% t(D2) %*% V2inv %*% R2))/(matTr(V1inv %*% R1) + matTr(V2inv %*% R2))
    hes2 = (matTr(V1inv %*% dV1 %*% V1inv %*% R1) + matTr(V2inv %*% dV2 %*% V2inv %*% R2))/(matTr(V1inv %*% R1) + matTr(V2inv %*% R2))
    val = - matTr(n1 * C1C1 + n2 *C2C2 ) + tmp * 2 * hes1 - tmp * hes2^2
    
    ##for REML
    if (lklMethod == "REML") {
      val = val +  matTr(C1C1) + matTr(C2C2)
    }
    return(val)
  }
  
  tau = newton(x0, lb, ub, f, g, h, tol=tol)$solution
  
  tmp = matTr(chol2inv(chol(tau * D1D1 + Ip)) %*% R1) + matTr(chol2inv(chol(tau * D2D2 + Ip)) %*% R2)
  
  se = ifelse((lklMethod == "REML"), (1/(N - 2*p)) * tmp, (1/N) * tmp)
  sg = tau * se
  
  return(list(s2e = se, s2g = sg, tau = tau))
}
