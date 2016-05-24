Taugoal <- function(x, k, data, scales, Tmax, controllo) {  
  ## x: beta,gamma
  beta <- x[1:k]
  gamma <- x[(k+1):length(x)]  
  dV <- dim(data$V)
  p <- dV[1]
  JL <- p*(p-1)/2
## SET STORAGE MODE OF y, x and V
  storage.mode(data$y) <- "double"
  storage.mode(data$x) <- "double"
  storage.mode(data$V) <- "double"
  resid <- vcrobresid(y=data$y, x=data$x, beta=beta)  
  n <- ncol(resid)
  Sigma <- Vprod(V=data$V, gamma=gamma) ## V0 is added automatically
  if (any(eigen(Sigma)$values <= 0))
    return(Tmax)
  RR <- rssr(resid=resid, Sigma=Sigma)
  scales <- doSsteppw(RR=RR, scale=scales, bb=controllo$bb, cc=controllo$tuning.chi, psi=controllo$psi, tol=controllo$rel.tol.scale, verbose=(controllo$trace.lev>2))
  if (controllo$psi=="optimal")
    scales <- scales*controllo$tuning.chi^2
  T <- doTausteppw(RR=RR, scale=scales, cc=controllo$tuning.psi, psi=controllo$psi)
  return(T)
}

eta0Taugoal <- function(beta, gamma, data, controllo) {
  v <- qchisq(seq(0.0001,0.9999,length=5000), 2)
  s0 <- doSstep(m=v, scale=1, bb=controllo$bb, cc=controllo$tuning.chi, psi=controllo$psi, tol=controllo$rel.tol.scale, verbose=(controllo$trace.lev>2))
  dV <- dim(data$V)
  p <- dV[1]
  JL <- p*(p-1)/2
## SET STORAGE MODE OF y, x and V
  storage.mode(data$y) <- "double"
  storage.mode(data$x) <- "double"
  storage.mode(data$V) <- "double"
  resid <- vcrobresid(y=data$y, x=data$x, beta=beta)  
  n <- ncol(resid)
  Sigma <- Vprod(V=data$V, gamma=gamma)
  RSR <- rsr(resid=resid, Sigma=Sigma)  
  eta0 <- doSstep(m=RSR/s0, scale=1, bb=controllo$bb, cc=controllo$tuning.chi, psi=controllo$psi, tol=controllo$rel.tol.scale, verbose=(controllo$trace.lev>2))
  return(eta0)
}

TaugoalDet <- function(x, k, data, scales, Tmax, controllo) {  
  ## x: beta,gamma
  beta <- x[1:k]
  gamma <- x[(k+1):length(x)]  
  dV <- dim(data$V)
  p <- dV[1]
  JL <- p*(p-1)/2
## SET STORAGE MODE OF y, x and V
  storage.mode(data$y) <- "double"
  storage.mode(data$x) <- "double"
  storage.mode(data$V) <- "double"
  resid <- vcrobresid(y=data$y, x=data$x, beta=beta)  
  n <- ncol(resid)
  Sigma <- Vprod(V=data$V, gamma=gamma) ## V0 is added automatically
  if (any(eigen(Sigma)$values <= 0))
    return(Tmax)
  detS <- sdet(Sigma)  
  RR <- rssr(resid=resid, Sigma=Sigma)
  scales <- doSsteppw(RR=RR, scale=scales, bb=controllo$bb, cc=controllo$tuning.chi, psi=controllo$psi, tol=controllo$rel.tol.scale, verbose=(controllo$trace.lev>2))
  if (controllo$psi=="optimal")
    scales <- scales*controllo$tuning.chi^2
  T <- doTausteppwDet(RR=RR, scale=scales, cc=controllo$tuning.psi, psi=controllo$psi, detS=detS)
  return(T)
}
