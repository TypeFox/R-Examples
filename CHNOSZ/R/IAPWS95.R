# functions for properties of water using
# the IAPWS-95 formulation (Wagner and Pruss, 2002)

IAPWS95.idealgas <- function(p, delta, tau) {
  ## the ideal gas part in the IAPWS-95 formulation
  # from Table 6.1 of Wagner and Pruss, 2002
  n <- c( -8.32044648201, 6.6832105268, 3.00632, 0.012436,
           0.97315, 1.27950, 0.96956, 0.24873 )
  gamma <- c( NA, NA, NA, 1.28728967, 
              3.53734222, 7.74073708, 9.24437796, 27.5075105 )
  # Equation 6.5
  phi <- function() log(delta) + n[1] + n[2]*tau + n[3]*log(tau) +
    sum( n[4:8] * log(1-exp(-gamma[4:8]*tau)) )
  # derivatives from Table 6.4
  phi.delta <- function() 1/delta+0+0+0+0
  phi.delta.delta <- function() -1/delta^2+0+0+0+0
  phi.tau <- function() 0+0+n[2]+n[3]/tau+sum(n[4:8]*gamma[4:8]*((1-exp(-gamma[4:8]*tau))^-1-1))
  phi.tau.tau <- function() 0+0+0-n[3]/tau^2-sum(n[4:8]*gamma[4:8]^2 * 
    exp(-gamma[4:8]*tau)*(1-exp(-gamma[4:8]*tau))^-2)
  phi.delta.tau <- function() 0+0+0+0+0
  return(get(p)())
}

IAPWS95.residual <- function(p, delta, tau) {
  ## the residual part in the IAPWS-95 formulation
  # from Table 6.2 of Wagner and Pruss, 2002
  c <- c(rep(NA,7),rep(1,15),rep(2,20),rep(3,4),4,rep(6,4),rep(NA,5))
  d <- c(1,1,1,2,2,3,4,1,1,1,2,2,3,4,
         4,5,7,9,10,11,13,15,1,2,2,2,3,4,
         4,4,5,6,6,7,9,9,9,9,9,10,10,12,
         3,4,4,5,14,3,6,6,6,3,3,3,NA,NA)
  t <- c(-0.5,0.875,1,0.5,0.75,0.375,1,4,6,12,1,5,4,2,
         13,9,3,4,11,4,13,1,7,1,9,10,10,3,
         7,10,10,6,10,10,1,2,3,4,8,6,9,8,
         16,22,23,23,10,50,44,46,50,0,1,4,NA,NA)
  n <- c( 0.12533547935523E-1, 0.78957634722828E1 ,-0.87803203303561E1 ,
          0.31802509345418   ,-0.26145533859358   ,-0.78199751687981E-2,
          0.88089493102134E-2,-0.66856572307965   , 0.20433810950965   ,
         -0.66212605039687E-4,-0.19232721156002   ,-0.25709043003438   ,
          0.16074868486251   ,-0.40092828925807E-1, 0.39343422603254E-6,
         -0.75941377088144E-5, 0.56250979351888E-3,-0.15608652257135E-4,
          0.11537996422951E-8, 0.36582165144204E-6,-0.13251180074668E-11,
         -0.62639586912454E-9,-0.10793600908932   , 0.17611491008752E-1,
          0.22132295167546   ,-0.40247669763528   , 0.58083399985759   ,
          0.49969146990806E-2,-0.31358700712549E-1,-0.74315929710341   ,
          0.47807329915480   , 0.20527940895948E-1,-0.13636435110343   ,
          0.14180634400617E-1, 0.83326504880713E-2,-0.29052336009585E-1,
          0.38615085574206E-1,-0.20393486513704E-1,-0.16554050063734E-2,
          0.19955571979541E-2, 0.15870308324157E-3,-0.16388568342530E-4,
          0.43613615723811E-1, 0.34994005463765E-1,-0.76788197844621E-1,
          0.22446277332006E-1,-0.62689710414685E-4,-0.55711118565645E-9,
         -0.19905718354408   , 0.31777497330738   ,-0.11841182425981   ,
         -0.31306260323435E2 , 0.31546140237781E2 ,-0.25213154341695E4 ,
         -0.14874640856724   , 0.31806110878444)
  alpha <- c(rep(NA,51),20,20,20,NA,NA)
  beta <- c(rep(NA,51),150,150,250,0.3,0.3)
  gamma <- c(rep(NA,51),1.21,1.21,1.25,NA,NA)
  epsilon <- c(rep(NA,51),1,1,1,NA,NA)
  a <- c(rep(NA,54),3.5,3.5)
  b <- c(rep(NA,54),0.85,0.95)
  B <- c(rep(NA,54),0.2,0.2)
  C <- c(rep(NA,54),28,32)
  D <- c(rep(NA,54),700,800)
  A <- c(rep(NA,54),0.32,0.32)
  # from Table 6.5
  i1 <- 1:7
  i2 <- 8:51
  i3 <- 52:54
  i4 <- 55:56
  # deriviatives of distance function
  Delta <- function(i) { Theta(i)^2 + B[i] * ((delta-1)^2)^a[i] }
  Theta <- function(i) { (1-tau) + A[i] * ((delta-1)^2)^(1/(2*beta[i])) }
  Psi <- function(i) { exp ( -C[i]*(delta-1)^2 - D[i]*(tau-1)^2 ) }
  dDelta.bi.ddelta <- function(i) { b[i]*Delta(i)^(b[i]-1)*dDelta.ddelta(i) }
  d2Delta.bi.ddelta2 <- function(i) { b[i]*( Delta(i)^(b[i]-1) * d2Delta.ddelta2(i) + 
    (b[i]-1)*Delta(i)^(b[i]-2)*dDelta.ddelta(i)^2 ) }
  dDelta.bi.dtau <- function(i) { -2*Theta(i)*b[i]*Delta(i)^(b[i]-1) }
  d2Delta.bi.dtau2 <- function(i) { 2*b[i]*Delta(i)^(b[i]-1) + 4*Theta(i)^2*b[i]*(b[i]-1)*Delta(i)^(b[i]-2) }
  d2Delta.bi.ddelta.dtau <- function(i) { -A[i]*b[i]*2/beta[i]*Delta(i)^(b[i]-1)*(delta-1) * 
    ((delta-1)^2)^(1/(2*beta[i])-1) - 2*Theta(i)*b[i]*(b[i]-1)*Delta(i)^(b[i]-2)*dDelta.ddelta(i)  }
  dDelta.ddelta <- function(i) { (delta-1) * ( A[i]*Theta(i)*2/beta[i]*((delta-1)^2)^(1/(2*beta[i])-1) +
    2*B[i]*a[i]*((delta-1)^2)^(a[i]-1) ) }
  d2Delta.ddelta2 <- function(i) { 1/(delta-1)*dDelta.ddelta(i) + (delta-1)^2 * (
    4*B[i]*a[i]*(a[i]-1)*((delta-1)^2)^(a[i]-2) + 2*A[i]^2*(1/beta[i])^2 *
      (((delta-1)^2)^(1/(2*B[i])-1))^2 + A[i]*Theta(i)*4/beta[i]*(1/(2*B[i])-1) *
        ((delta-1)^2)^(1/(2*beta[i])-2) ) }
  # derivatives of exponential function
  dPsi.ddelta <- function(i) { -2*C[i]*(delta-1)*Psi(i) }
  d2Psi.ddelta2 <- function(i) { ( 2*C[i]*(delta-1)^2 - 1 ) * 2*C[i]*Psi(i) }
  dPsi.dtau <- function(i) { -2*D[i]*(tau-1)*Psi(i) }
  d2Psi.dtau2 <- function(i) { (2*D[i]*(tau-1)^2 - 1) * 2*D[i]*Psi(i) }
  d2Psi.ddelta.dtau <- function(i) { 4*C[i]*D[i]*(delta-1)*(tau-1)*Psi(i) }
  # dimensionless Helmholtz free energy and derivatives
  phi <- function() {
    sum(n[i1]*delta^d[i1]*tau^t[i1]) +
    sum(n[i2]*delta^d[i2]*tau^t[i2]*exp(-delta^c[i2])) +
    sum(n[i3]*delta^d[i3]*tau^t[i3] *
      exp( -alpha[i3]*(delta-epsilon[i3])^2 - beta[i3]*(tau-gamma[i3])^2 ) ) +
    sum(n[i4]*Delta(i4)^b[i4]*delta*Psi(i4))
  }
  phi.delta <- function() {
    sum(n[i1]*d[i1]*delta^(d[i1]-1)*tau^t[i1]) +
    sum(n[i2]*exp(-delta^c[i2])*(delta^(d[i2]-1)*tau^t[i2]*(d[i2]-c[i2]*delta^c[i2]))) +
    sum(n[i3]*delta^d[i3]*tau^t[i3] *
      exp( -alpha[i3]*(delta-epsilon[i3])^2 - beta[i3]*(tau-gamma[i3])^2 ) * 
        (d[i3]/delta - 2 * alpha[i3]*(delta-epsilon[i3])) ) +
    sum(n[i4] * ( Delta(i4)^b[i4] * (Psi(i4)+delta*dPsi.ddelta(i4)) + dDelta.bi.ddelta(i4)*delta*Psi(i4) ) ) 
  }
  phi.delta.delta <- function() {
    sum(n[i1]*d[i1]*(d[i1]-1)*delta^(d[i1]-2)*tau^t[i1]) +
    sum(n[i2]*exp(-delta^c[i2])*(delta^(d[i2]-2)*tau^t[i2]*((d[i2]-c[i2]*delta^c[i2]) * 
      (d[i2]-1-c[i2]*delta^c[i2])-c[i2]^2*delta^c[i2]))) +
    sum(n[i3]*tau^t[i3]*exp(-alpha[i3]*(delta-epsilon[i3])^2 - beta[i3]*(tau-gamma[i3])^2) * (
      -2*alpha[i3]*delta^d[i3]+4*alpha[i3]^2*delta^d[i3]*(delta-epsilon[i3])^2 -
       4*d[i3]*alpha[i3]*delta^(d[i3]-1)*(delta-epsilon[i3])+d[i3]*(d[i3]-1)*delta^(d[i3]-2) ) ) +
    sum(n[i4]*( Delta(i4)^b[i4]*(2*dPsi.ddelta(i4)+delta*d2Psi.ddelta2(i4)) + 
      2*dDelta.bi.ddelta(i4)*(Psi(i4)+delta*dPsi.ddelta(i4)) + d2Delta.bi.ddelta2(i4)*delta*Psi(i4) ) )
  }
  phi.tau <- function() {
    sum(n[i1]*t[i1]*delta^d[i1]*tau^(t[i1]-1)) +
    sum(n[i2]*t[i2]*delta^d[i2]*tau^(t[i2]-1)*exp(-delta^c[i2])) +
    sum(n[i3]*delta^d[i3]*tau^t[i3]*exp(-alpha[i3]*(delta-epsilon[i3])^2-beta[i3]*(tau-gamma[i3])^2) * 
      (t[i3]/tau-2*beta[i3]*(tau-gamma[i3]))) +
    sum(n[i4]*delta*(dDelta.bi.dtau(i4)*Psi(i4)+Delta(i4)^b[i4]*dPsi.dtau(i4)))
  }
  phi.tau.tau <- function() {
    sum(n[i1]*t[i1]*(t[i1]-1)*delta^d[i1]*tau^(t[i1]-2)) +
    sum(n[i2]*t[i2]*(t[i2]-1)*delta^d[i2]*tau^(t[i2]-2)*exp(-delta^c[i2])) +
    sum(n[i3]*delta^d[i3]*tau^t[i3]*exp(-alpha[i3]*(delta-epsilon[i3])^2-beta[i3]*(tau-gamma[i3])^2) * 
      (((t[i3]/tau)-2*beta[i3]*(tau-gamma[i3]))^2-t[i3]/tau^2-2*beta[i3])) +
    sum(n[i4]*delta*(d2Delta.bi.dtau2(i4)*Psi(i4)+2*dDelta.bi.dtau(i4)*dPsi.dtau(i4) + 
      Delta(i4)^b[i4]*d2Psi.dtau2(i4)))
  }
  phi.delta.tau <- function() {
    sum(n[i1]*d[i1]*t[i1]*delta^(d[i1]-1)*tau^(t[i1]-1)) +
    sum(n[i2]*t[i2]*delta^(d[i2]-1)*tau^(t[i2]-1)*(d[i2]-c[i2]*delta^c[i2])*exp(-delta^c[i2])) +
    sum(n[i3]*delta^d[i3]*tau^t[i3]*exp(-alpha[i3]*(delta-epsilon[i3])^2-beta[i3]*(tau-gamma[i3])^2) * 
      ((d[i3]/delta)-2*alpha[i3]*(delta-epsilon[i3]))*(t[i3]/tau-2*beta[i3]*(tau-gamma[i3])) ) +
    sum(n[i4]*(Delta(i4)^b[i4]*(dPsi.dtau(i4)+delta*d2Psi.ddelta.dtau(i4)) + 
      delta*dDelta.bi.ddelta(i4)*dPsi.dtau(i4)+dDelta.bi.dtau(i4) * (Psi(i4)+delta*dPsi.ddelta(i4)) +
      d2Delta.bi.ddelta.dtau(i4)*delta*Psi(i4) )) 
  }
  return(get(p)())
}

IAPWS95 <- function(property,T=298.15,rho=1000) {
  ## the IAPWS-95 formulation for ordinary water substance
  ## Wagner and Pruss, 2002
  property <- tolower(property)
  # triple point
  T.triple <- 273.16 # K
  P.triple <- 611.657 # Pa
  rho.triple.liquid <- 999.793
  rho.triple.vapor <- 0.00485458
  # normal boiling point
  T.boiling <- 373.124
  P.boiling <- 0.101325
  rho.boiling.liquid <- 958.367
  rho.boiling.vapor <- 0.597657
  # critical point constants
  T.critical <- 647.096 # K
  rho.critical <- 322 # kg m-3
  # specific and molar gas constants
  R <- 0.46151805 # kJ kg-1 K-1
  # R.M <- 8.314472 # J mol-1 K-1
  # molar mass
  M <- 18.015268 # g mol-1
  ## define functions idealgas and residual, supplying arguments delta and tau
  idealgas <- function(p) IAPWS95.idealgas(p, delta, tau)
  residual <- function(p) IAPWS95.residual(p, delta, tau)
  ## relation of thermodynamic properties to Helmholtz free energy
  a <- function() {
    x <- idealgas('phi')+residual('phi')
    return(x*R*T)
  }
  # Table 6.3 
  p <- function() {
    x <- 1 + delta*residual('phi.delta')
    return(x*rho*R*T/1000)  # for MPa
  }
  s <- function() {
    x <- tau * (idealgas('phi.tau')+residual('phi.tau'))-idealgas('phi')-residual('phi')
    return(x*R)
  }
  u <- function() {
    x <- tau * (idealgas('phi.tau')+residual('phi.tau'))
    return(x*R*T)
  }
  h <- function() {
    x <- 1 + tau * (idealgas('phi.tau')+residual('phi.tau')) + delta*residual('phi.delta')
    return(x*R*T)
  }
  g <- function() {
    x <- 1 + idealgas('phi') + residual('phi') + delta*residual('phi.delta')
    return(x*R*T)
  }
  cv <- function() {
    x <- -tau^2*(idealgas('phi.tau.tau')+residual('phi.tau.tau'))
    return(x*R)
  }
  cp <- function() {
    x <- -tau^2*(idealgas('phi.tau.tau')+residual('phi.tau.tau')) +
         (1+delta*residual('phi.delta')-delta*tau*residual('phi.delta.tau'))^2 /
         (1+2*delta*residual('phi.delta')+delta^2*residual('phi.delta.delta'))
    return(x*R)
  }
# 20090420 speed of sound calculation is incomplete
# (delta.liquid and drhos.dT not visible)
#  cs <- function() {
#    x <- -tau^2*(idealgas('phi.tau.tau')+residual('phi.tau.tau')) +
#         (1+delta*residual('phi.delta')-delta*tau*residual('phi.delta.tau'))^2 /
#         (1+2*delta*residual('phi.delta')+delta^2*residual('phi.delta.delta')) *
#         ((1+delta.liquid*residual('phi.delta')-delta.liquid*tau*residual('phi.tau.tau'))-rho.critical/(R*delta.liquid)*drhos.dT)
#    return(x*R)
#  }
  w <- function() {
    x <- 1 + 2*delta*residual('phi.delta') + delta^2*residual('phi.delta.delta') - 
         (1+delta*residual('phi.delta')-delta*tau*residual('phi.delta.tau'))^2 /
         tau^2*(idealgas('phi.tau.tau')+residual('phi.tau.tau'))
    return(sqrt(x*R*T))
  }
  mu <- function() {
    x <- -(delta*residual('phi.delta')+delta^2*residual('phi.delta.delta')+delta*tau*residual('phi.delta.tau')) /
          ( ( 1+delta*residual('phi.delta')-delta*tau*residual('phi.delta.tau')^2 ) - tau^2 *
          (idealgas('phi.tau.tau')+residual('phi.tau.tau'))*(1+2*delta*residual('phi.delta')+delta^2*residual('phi.delta.delta')) ) 
    return(x/(R*rho))
  }
  ## run the calculations
  ww <- NULL
  my.T <- T
  my.rho <- rho
  for(j in 1:length(property)) {
    t <- numeric()
    for(i in 1:length(my.T)) {
      T <- my.T[i]
      rho <- my.rho[i]
      # Equation 6.4
      delta <- rho / rho.critical
      tau <- T.critical / T
      t <- c(t,get(property[j])())
    }
    t <- data.frame(t)
    if(j==1) ww <- t else ww <- cbind(ww,t)
  }
  colnames(ww) <- property
  return(ww)
}

