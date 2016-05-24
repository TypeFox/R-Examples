approxVarEst <-
function(se0, sg0, D1, D2, r, n1, n2, control=NULL){
  if (is.null(control)) {
    tol = 0.01
    s2profile = "se"
    lklMethod = "REML"
  } else {
    tol = control$tol
    s2profile = control$s2profile
    lklMethod = control$lklMethod
  }
  
  matTr <- function(z) sum(diag(z))
  
  p = dim(D1)[1]
  Ip = diag(rep(1, p))
  
  n = n1 + n2
  N = n * p
  
  ##residual matrices\n
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
  
  gap = 1
  sg = sg0
  se = se0
  cnt = 0
  ## Whether it's ML or REML
  lklConst = ifelse(lklMethod == "REML", N - 2*p, N)
  
  while (gap > tol) {
    sg0 = sg
    se0 = se
    cnt = cnt + 1
    
    tmp0 = matTr(solve(sg * D1D1 + se * Ip) %*% R1) + matTr(solve(sg * D2D2 + se * Ip) %*% R2)
    
    if (s2profile == "se") {     
      se =  tmp0 / lklConst
      
      tmp1 = cov(t(r[, 1:n1]))
      tmp = min(diag(tmp1))
      tmp = min(se, tmp)
      tmp1 = tmp1 - diag(rep(tmp, p))
      tmp1 = solve(D1) %*% tmp1
      
      tmp2 = cov(t(r[, (n1 + 1):n]))
      tmp = min(diag(tmp2))
      tmp = min(se, tmp)
      tmp2 = tmp2 - diag(rep(tmp, p))
      tmp2 = solve(D2) %*% tmp2
      
      sg = (mean(diag(tmp1)) + mean(diag(tmp2)))/2
      tau = sg/se     
    } else {     
      sg = tmp0/lklConst
      
      tmp1 = cov(t(r[, 1:n1]))
      tmp = min(diag(tmp1))
      tmp = min(sg, tmp)
      tmp1 = tmp1 - diag(rep(tmp, p))
      
      tmp2 = cov(t(r[, (n1 + 1):n]))
      tmp = min(diag(tmp2))
      tmp = min(sg, tmp)
      tmp2 = tmp2 - diag(rep(tmp, p))
      
      se = (mean(diag(tmp1)) + mean(diag(tmp2)))/2
      tau = se/sg
    }
    
    gap = abs(sg - sg0) + abs(se - se0)
  }
  
  return(list(s2e = se, s2g = sg, tau=tau, finalgap = gap, niter = cnt))
}
