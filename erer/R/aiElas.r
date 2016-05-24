aiElas <- function(z, ...) {
  if (!inherits(z, "aiFit")) {
      stop("\nPlease provide an object from 'aiStaFit' or 'aiDynFit'.\n") }
  est <- z$est;      
  cof <- coef(est);  vco   <- vcov(est);  df <- df.residual(est)	
  nE  <- z$nExoge;   nS    <- z$nShare;   nP <- z$nParam; nO <- z$nOmit 
  if (inherits(z, "aiStaFit")) {e.share <- z$share  }
  if (inherits(z, "aiDynFit")) {e.share <- z$w$share}
  av <- colMeans(z$y[, e.share])
  av.sh <- av[c(e.share[-nO], e.share[nO])]  
 
  # coefficient for the beta and gamma in the omitted equation
  bt.co <- bt.vc <- gm.co <- gm.vc <- 0
  for ( i in 1:(nS-1)) { 
    bt.co <- bt.co - cof[nE + nP*(i-1)]
    gm.co <- gm.co - cof[nP*i] 
    for ( j in 1:(nS-1)) {
      bt.vc <- bt.vc + vco[nE + nP*(i-1), nE + nP*(j-1)]
      gm.vc <- gm.vc + vco[nP*i, nP*j]       
    }
  }
  
  # ex = Expenditure elasticity, including that for omitted equation
  e.ex <- v.ex <- NULL
  for(i in 1:(nS-1)){
      e.ex[i] <- 1 + cof[nE + nP*(i-1)] / av.sh[i]
      v.ex[i] <- vco[nE + nP*(i-1), nE + nP*(i-1)] / (av.sh[i]^2)
  }
  e.ex[nS] <- 1 + bt.co / av.sh[nS]
  v.ex[nS] <- bt.vc / (av.sh[nS]^2)

  # Marshallian and Hicksian elasticiity: e = estimate, v = variance
  e.ma <- v.ma <- e.hi <- v.hi <-matrix(NA, nS, nS)
  for (i in 1:(nS-1)) {
    for (j in 1:nS) {
      delta <- ifelse(i==j, 1, 0)
      e.ma[i,j] <- (-delta + cof[(i-1)*nP+nE+j] / av.sh[i] -
                     cof[nE + nP*(i-1)] * av.sh[j]/av.sh[i])
      e.hi[i,j] <- -delta + cof[(i-1)*nP+nE+j] / av.sh[i] + av.sh[j]
      v.ma[i,j] <- ( vco[(i-1)*nP+nE+j, (i-1)*nP+nE+j] / (av.sh[i]^2) +
                     vco[(i-1)*nP+nE,   (i-1)*nP+nE  ] * (av.sh[j]^2) / (av.sh[i]^2) - 
                     vco[(i-1)*nP+nE+j, (i-1)*nP+nE  ] * av.sh[j] * 2 / (av.sh[i]^2) )
      v.hi[i,j] <- vco[(i-1)*nP+nE+j,(i-1)*nP+nE+j] / (av.sh[i]^2)
    }
  }
  for (j in 1:(nS-1)) {
      e.ma[nS, j] <- cof[j*nP] / av.sh[nS] - av.sh[j] * bt.co / av.sh[nS]
      e.hi[nS, j] <- cof[j*nP] / av.sh[nS] + av.sh[j]
      v.ma[nS, j] <- (vco[j*nP, j*nP] + bt.vc *(av.sh[j]^2)) / (av.sh[nS]^2)
      v.hi[nS, j] <- vco[j*nP, j*nP] / (av.sh[nS]^2)
  }
  e.ma[nS, nS] <-  -1 + gm.co /av.sh[nS] - bt.co    # omitted equation
  e.hi[nS, nS] <-  -1 + gm.co /av.sh[nS] + av.sh[nS]
  v.ma[nS, nS] <-  gm.vc / (av.sh[nS]^2) + bt.vc
  v.hi[nS, nS] <-  gm.vc / (av.sh[nS]^2)

  # Expenditure: t-ratio, p value, all
  t.ex <- e.ex/sqrt(v.ex)
  t.ma <- e.ma/sqrt(v.ma)
  t.hi <- e.hi/sqrt(v.hi) 
  p.ex <- 2*( 1- pt(abs(t.ex),df) )
  p.ma <- 2*( 1- pt(abs(t.ma),df) )
  p.hi <- 2*( 1- pt(abs(t.hi),df) )    
  ex <- data.frame(cbind(e.ex, sqrt(v.ex), t.ex, p.ex)) # expenditure
  rownames(ex) <- names(av.sh)
  expen <- bsTab(ex, ...)
  colnames(expen) <- c("Elas.expen", "Estimate")
 
 # Marshallian and Hicksian: t, p, all
  mars <- hick <- data.frame(matrix(NA, nrow=2*nS, ncol=nS))
  for (i in 1:nS) {  
    ma <- data.frame(cbind(e.ma[,i], sqrt(v.ma)[,i], t.ma[,i], p.ma[,i]))
    hi <- data.frame(cbind(e.hi[,i], sqrt(v.hi)[,i], t.hi[,i], p.hi[,i])) 
    rownames(ma) <- rownames(hi) <- names(av.sh)       
    mars[, i] <- bsTab(ma, ...)[,2]
    hick[, i] <- bsTab(hi, ...)[,2]     
  }
  colnames(mars) <- colnames(hick) <- names(av.sh)
  marsh <- cbind(Elas.Marsh = bsTab(ma, ...)[1], mars)
  hicks <- cbind(Elas.Hicks = bsTab(hi, ...)[1], hick)    
  
  result <- listn(name = names(av.sh), expen, marsh, hicks)
  return(result)
}  