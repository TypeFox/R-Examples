
ram_calc = function(par,SampCov22,A,S,F,SampMean){

  #mult = extractMatrices(fit.growth)
  #A <- mult$A
 # S <- mult$S
 # F <- mult$F


  ret <- list()

  #coefs = coef(fit.growth)
 # par = c(coefs[1:2],coefs[10:11],coefs[3:9])

  A2 <- A
  S2 <- S
  # doesn't work for
  for(i in 1:length(par)){
    A2[A2== i] <- par[i]
    S2[S2== i] <- par[i]
  }

  #A2[A.fixed] <- A.est[A.fixed]
  #S2[S.fixed] <- S.est[S.fixed]




  #ImpCov = F %*% solve(I-A2) %*% S2 %*% t(solve(I-A2)) %*% t(F)


  nncol = which(colnames(A2) == "1")
  m = A2[-nncol,"1"]
  m.pars = A[-nncol,"1"]
  #m.pars = A[-nncol,"1"]
  A.pars = A[-nncol,-nncol]
  A2 = A2[-nncol,-nncol]
  F= F[-nncol,-nncol]
  dd = SampMean
  S2 = S2[-nncol,-nncol]
  S.pars = S[-nncol,-nncol]
  I = diag(nrow(A2))
  mu = F %*% solve(I - A2) %*% m



  #I = diag(nrow(Areg))
  ImpCov = F %*% solve(I-A2) %*% S2 %*% t(solve(I-A2)) %*% t(F)

  #grad.out <- rep(0,length(par))

  B = solve(I - A2)
  C = diag(nrow(ImpCov)) - (solve(ImpCov) %*% SampCov22)
  E = B %*% S2 %*% t(B)
  b = dd - mu

  lik = log(((2*pi)^length(b)) * det(ImpCov)) + trace(solve(ImpCov) %*% SampCov22) +  t(b) %*% solve(ImpCov) %*% b
  #lik= log(det(ImpCov)) + trace(SampCov22 %*% solve(ImpCov)) - log(det(SampCov22))  - 4


  #A.iter <- max(A)



  ret$lik <- lik
  ret$ImpCov <- ImpCov
  ret$S2 <- S2
  ret$A2 <- A2
  ret$m <- m
  ret$m.pars <- m.pars
  ret$A.pars <- A.pars
  ret$S.pars <- S.pars
  ret$F <- F
  ret$mu <- mu
  #lik;ImpCov;S2;A2
  ret
}
