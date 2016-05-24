profile.newton.sg <-
function(x0, D1, D2, r, n1, n2, control=NULL) {
  if(is.null(control)){   
    lklMethod = 'REML'  
    lb = 0.01         	
    ub = 100 		       	
    s2profile = 'sg' 		
    tol = 0.01          
  } else{    
    lklMethod = control$lklMethod
    lb = control$lb
    ub = control$ub
    s2profile = control$s2profile    
    tol = control$tol
  }  
  
  p = nrow(D1)
  Ip = diag(rep(1, p), p, p)
  n = n1 + n2
  
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
  
  ##matrix trace function
  matTr <- function(z) sum(diag(z))
  
  #obj fn / the same as before
  f <- function(tau) {
    V1 = D1D1 + tau * Ip
    V2 = D2D2 + tau * Ip
    V1inv = chol2inv(chol(V1))
    V2inv = chol2inv(chol(V2))
    
    tmp = ifelse((lklMethod == "REML"), N - 2* p, N)
    
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
  #-----
  
  #-----
  #gradient fn
  
  g <- function(tau) {
    V1 = D1D1 + tau * Ip
    V2 = D2D2 + tau * Ip
    V1inv = chol2inv(chol(V1))
    V2inv = chol2inv(chol(V2))
    C1 = V1inv %*% V1inv
    C2 = V2inv %*% V2inv
    
    tmp = ifelse((lklMethod == "REML"), N - 2* p, N)
    
    val = n1 * matTr(V1inv) + n2 * matTr(V2inv) - tmp * (matTr(C1 %*% R1) + matTr(C2 %*% R2))/(matTr(V1inv %*% R1) + matTr(V2inv %*% R2))
    ##for REML
    if (lklMethod == "REML") {
      HC = t(D1) %*% V1inv %*% D1 
      HT = t(D2) %*% V2inv %*% D2
      val = val - matTr(chol2inv(chol(HC)) %*% t(D1) %*% C1 %*% D1 ) - matTr( chol2inv(chol(HT))%*% t(D2) %*% C2 %*% D2)
    }
    
    return(val)
  }
  #-----
  
  #-----
  #hessian fn
  
  h <- function(tau) {
    V1 = D1D1 + tau * Ip
    V2 = D2D2 + tau * Ip
    V1inv = chol2inv(chol(V1))
    V2inv = chol2inv(chol(V2))
    C1 = V1inv %*% V1inv
    C2 = V2inv %*% V2inv
    
    tmp = ifelse((lklMethod == "REML"), N - 2*p, N)
    
    hes1 = (matTr(C1 %*% V1inv %*% R1) + matTr(C2 %*% V2inv %*% R2))/(matTr(V1inv %*% R1) + matTr(V2inv %*% R2))
    hes2 = (matTr(C1 %*% R1) + matTr(C2 %*% R2))/(matTr(V1inv %*% R1) + matTr(V2inv %*% R2))
    val = -matTr(n1 * C1 + n2 * C2) + tmp * 2 * hes1 - tmp * hes2^2
    
    ##for REML
    if (lklMethod == "REML") {
      HC = t(D1) %*% V1inv %*% D1 
      HT = t(D2) %*% V2inv %*% D2
      HCinv = chol2inv(chol(HC)) 
      HTinv = chol2inv(chol(HT))
      val = val - matTr( HCinv %*% t(D1) %*% C1 %*% D1 ) - matTr( HTinv %*% t(D2) %*% C2 %*% D2)
      
      val = val + 2 * matTr(HCinv %*% t(D1) %*% C1 %*% V1inv %*% D1 )  + 2*matTr(HTinv %*% t(D2) %*% C2 %*% V2inv %*% D2)
    }
    return(val)
  }
  
  tau = newton(x0, lb, ub, f, g, h, tol=tol)$solution  
  
  tmp = matTr(chol2inv(chol(D1D1 + tau * Ip)) %*% R1) + matTr(chol2inv(chol(D2D2 + tau * Ip)) %*% R2)
  
  sg = ifelse((lklMethod == "REML"), (1/(N - 2*p)) * tmp, (1/N) * tmp)
  se = tau * sg
  
  return(list(s2e = se, s2g = sg, tau = tau))
}
