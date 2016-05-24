# CHNOSZ/util.water.R

WP02.auxiliary <- function(property='rho.liquid',T=298.15) {
  # auxiliary equations for liquid-vapor phase boundary
  # from Wagner and Pruss, 2002
  # critical point
  T.critical <- 647.096 # K
  P.critical <- 22.064 # MPa
  rho.critical <- 322 # kg m-3

  if(property %in% c("P.sigma","dP.sigma.dT")) {
    # vapor pressure
    V <- 1 - T / T.critical # theta (dimensionless)
    a1 <- -7.85951783
    a2 <- 1.84408259
    a3 <- -11.7866497
    a4 <- 22.6807411
    a5 <- -15.9618719
    a6 <- 1.80122502
    ln.P.sigma.P.critical <- T.critical / T *
      ( a1*V + a2*V^1.5 + a3*V^3 + a4*V^3.5 + a5*V^4 + a6*V^7.5 ) 
    P.sigma <- P.critical * exp(ln.P.sigma.P.critical)
    if(property=="dP.sigma.dT") out <- - P.sigma / T * ( ln.P.sigma.P.critical +
      a1 + 1.5*a2*V^0.5 + 3*a3*V^2 + 3.5*a4*V^2.5 + 4*a5*V^3 + 7.5*a6*V^6.5 )
    else out <- P.sigma
  } else if(property=="rho.liquid") {
    # saturated liquid density
    V <- 1 - T / T.critical
    b1 <- 1.99274064
    b2 <- 1.09965342
    b3 <- -0.510839303
    b4 <- -1.75493479
    b5 <- -45.5170352
    b6 <- -6.74694450E5
    rho.liquid <- rho.critical * ( 
      1 + b1*V^(1/3) + b2*V^(2/3) + b3*V^(5/3) + b4*V^(16/3) + b5*V^(43/3) + b6*V^(110/3) )
    out <- rho.liquid
  } else if(property=="rho.vapor") {
  # saturated vapor density
    V <- 1 - T / T.critical
    c1 <- -2.03150240
    c2 <- -2.68302940
    c3 <- -5.38626492
    c4 <- -17.2991605
    c5 <- -44.7586581
    c6 <- -63.9201063
    rho.vapor <- rho.critical * exp (
      c1*V^(2/6) + c2*V^(4/6) + c3*V^(8/6) + c4*V^(18/6) + c5*V^(37/6) + c6*V^(71/6) )
    out <- rho.vapor
  } else stop(paste('i can not calculate',property))
  return(out)
}

# return a density in kg m-3
# corresponding to the given pressure (bar) and temperature (K)
rho.IAPWS95 <- function(T=298.15, P=1, state="", trace=0) {
  # function for which to find a zero
  dP <- function(rho, T, P.MPa) IAPWS95("P", rho=rho, T=T)[, 1] - P.MPa
  # convert bar to MPa
  P.MPa <- convert(P, "MPa")
  rho <- numeric()
  T.critical <- 647.096 # K
  P.critical <- 22.064 # MPa
  for(i in 1:length(T)) {
    Psat <- WP02.auxiliary("P.sigma", T[i])
    if(T[i] > T.critical) {
      # above critical temperature
      interval <- c(0.1, 1)
      extendInt <- "upX"
      if(trace > 0) msgout("supercritical (T) ")
    } else if(P.MPa[i] > P.critical) {
      # above critical pressure
      rho.sat <- WP02.auxiliary("rho.liquid", T=T[i])
      interval <- c(rho.sat, rho.sat + 1)
      extendInt <- "upX"
      if(trace > 0) msgout("supercritical (P) ")
    } else if(P.MPa[i] <= 0.9999*Psat) {
      # steam
      rho.sat <- WP02.auxiliary("rho.vapor", T=T[i])
      interval <- c(rho.sat*0.1, rho.sat)
      extendInt <- "upX"
      if(trace > 0) msgout("steam ")
    } else if(P.MPa[i] >= 1.00005*Psat) {
      # water
      rho.sat <- WP02.auxiliary("rho.liquid", T=T[i])
      interval <- c(rho.sat, rho.sat + 1)
      extendInt <- "upX"
      if(trace > 0) msgout("water ")
    } else if(!state %in% c("liquid", "vapor")) {
      # we're close to the saturation curve;
      # calculate rho and G for liquid and vapor and return rho for stable phase
      if(trace > 0) msgout("close to saturation; trying liquid and vapor\n")
      rho.liquid <- rho.IAPWS95(T[i], P[i], state="liquid", trace=trace)
      rho.vapor <- rho.IAPWS95(T[i], P[i], state="vapor", trace=trace)
      G.liquid <- IAPWS95("G", rho=rho.liquid, T=T[i])
      G.vapor <- IAPWS95("G", rho=rho.vapor, T=T[i])
      if(G.liquid < G.vapor) {
        this.rho <- rho.liquid 
        if(trace > 0) msgout(paste0("G.liquid(", G.liquid, ") < G.vapor(", G.vapor, ")\n"))
      } else {
        this.rho <- rho.vapor
        if(trace > 0) msgout(paste0("G.vapor(", G.vapor, ") < G.liquid (", G.liquid, ")\n"))
      }
      rho <- c(rho, this.rho)
      next
    } else {
      # we are looking at a specific state
      if(trace > 0) msgout(paste("specified state:", state, " "))
      if(state=="vapor") rho0 <- WP02.auxiliary("rho.vapor", T[i])
      else if(state=="liquid") rho0 <- WP02.auxiliary("rho.liquid", T[i])
      # a too-big range may cause problems e.g.
      # interval <- c(rho0*0.9, rho0*1.1) fails for T=253.15, P=1
      interval <- c(rho0*0.95, rho0*1.05)
      # if P on the initial interval are both higher or lower than target P,
      # set the direction of interval extension
      P.init <- IAPWS95("P", rho=interval, T=c(T[i], T[i]))[, 1]
      if(all(P.init < P.MPa[i])) extendInt <- "downX"
      else if(all(P.init > P.MPa[i])) extendInt <- "upX"
      else extendInt <- "yes"
    }
    if(trace > 0) msgout(paste0("T=", T[i], " P=", P[i], " rho=[", interval[1], ",", interval[2], "]\n"))
    this.rho <- try(uniroot(dP, interval, extendInt=extendInt, trace=trace, T=T[i], P.MPa=P.MPa[i])$root, TRUE)
    if(!is.numeric(this.rho)) {
      warning("rho.IAPWS95: problems finding density at ", T[i], " K and ", P[i], " bar", call.=FALSE)
      this.rho <- NA
    }
    rho <- c(rho, this.rho)
  }
  return(rho)
}

water.AW90 <- function(T=298.15,rho=1000,P=0.1) {
  # Equations for the dielectric constant of water
  # from Archer and Wang, 1990
  # T in K
  # rho in kg m-3
  # p in MPa

  # Table 2
  b <- c(-4.044525E-2, 103.6180   , 75.32165   ,
         -23.23778   ,-3.548184   ,-1246.311   ,
         263307.7    ,-6.928953E-1,-204.4473)
  alpha <- 18.1458392E-30 # m^3
  #alpha <- 14.7E-30
  mu <- 6.1375776E-30 # C m
  N.A <- 6.0221367E23 # mol-1
  k <- 1.380658E-23 # Boltzmann constant, J K-1
  M <- 0.0180153 # kg mol-1
  rho.0 <- 1000 # kg m-3
  # Equation 1
  epsilon.0 <- 8.8541878E-12 # permittivity of vacuum, C^2 J-1 m-1
  #epsfun.lhs <- function(e) (e-1)*(2*e+1)/(9*e)
  epsfun.rhs <- function(T,V.m) N.A*(alpha+mufun()/(3*epsilon.0*k*T))/(3*V.m)
  #epsfun <- function(e,T,V.m) epsfun.lhs(e) - epsfun.rhs(T,V.m)
  mufun <- function() gfun()*mu^2
  gfun <- function() rhofun()*rho/rho.0 + 1
  # Equation 3
  rhofun <- function() b[1]*P*T^-1 + b[2]*T^-0.5 + b[3]*(T-215)^-1 +
    b[4]*(T-215)^-0.5 + b[5]*(T-215)^-0.25 +
    exp(b[6]*T^-1 + b[7]*T^-2 + b[8]*P*T^-1 + b[9]*P*T^-2)
  epsilon <- function(T,rho) {
    #tu <- try(uniroot(epsfun,c(1E-1,1E3),T=T,V.m=M/rho)$root,TRUE)
    epspoly <- function() epsfun.rhs(T,V.m=M/rho)
    tu <- (9*epspoly() + 1 + ((9*epspoly()+1)*(9*epspoly()+1) + 8)^0.5) / 4 #Marc Neveu added 4/24/2013
    if(!is.numeric(tu)) {
      warning('water.AW90: no root for density at ',T,' K and ',rho,' kg m-3.',call.=FALSE,immediate.=TRUE)
      tu <- NA
    }
    return(tu)
  }
  # get things the right length
  our.T <- T; our.rho <- rho; our.P <- P
  t <- numeric()
  for(i in 1:length(our.T)) {
    T <- our.T[i]
    rho <- our.rho[i]
    P <- our.P[i]
    t <- c(t,epsilon(T,rho))
  }
  return(t)
}
