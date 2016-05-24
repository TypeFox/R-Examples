# CHNOSZ/water.R
# calculate thermodynamic and electrostatic properties of H2O
# 20061016 jmd

water <- function(property = NULL, T = get("thermo")$opt$Tr, P = "Psat") {
  # calculate the properties of liquid H2O as a function of T and P
  # T in Kelvin, P in bar
  if(is.null(property)) stop('property was NULL')
  # make T and P equal length
  if(!identical(P, "Psat")) {
    if(length(P) < length(T)) P <- rep(P, length.out = length(T))
    else if(length(T) < length(P)) T <- rep(T, length.out = length(P))
  }
  # turn 273.15 K to 273.16 K (needed for water.SUPCRT92 at Psat)
  T[T == 273.15] <- 273.16
  # this tells us to do the calculations using code taken from SUPCRT
  do.supcrt <- get("thermo")$opt$water != "IAPWS95"
  if(do.supcrt) {
    # get the values of the properties using SUPCRT92
    w.out <- water.SUPCRT92(property, T, P)
    return(w.out)
  } else {
    # here we get properties using IAPWS-95 
    w.out <- water.IAPWS95(property, T, P)
    # normalize the names to use upper case (expected by subcrt())
    iprop <- match(tolower(property), tolower(water.props("IAPWS95")))
    colnames(w.out) <- water.props("IAPWS95")[iprop]
    return(w.out)
  }
}

water.props <- function(formulation=get("thermo")$opt$water) {
  # return the names of properties that are available in SUPCRT92 or IAPWS95
  # added 20130212 jmd
  if(formulation=="SUPCRT92")
    props <- c("A", "G", "S", "U", "H", "Cv", "Cp",
    "Speed", "alpha", "beta", "diel", "visc",
    "tcond", "surten", "tdiff", "Prndtl", "visck", "albe",
    "ZBorn", "YBorn", "QBorn", "daldT", "XBorn",
    "V", "rho", "Psat", "E", "kT")
  else if(formulation=="IAPWS95")
    props <- c("A", "G", "S", "U", "H", "Cv", "Cp",
    "Speed", "diel",
    "YBorn", "QBorn", "XBorn", "NBorn", "UBorn",
    "V", "rho", "Psat", "de.dT", "de.dP", "P")
  return(props)
}

water.SUPCRT92 <- function(property, T=298.15, P=1) {
  ### interface to H2O92D.f : FORTRAN subroutine taken from 
  ### SUPCRT92 for calculating the thermodynamic and 
  ### electrostatic properties of H2O. 
  ## we restrict the calculations to liquid water
  ## except for getting Psat (vapor-liquid saturation 
  ## pressure as a function of T>100 C). 20071213 jmd
  # check for availability of properties
  iprop <- match(tolower(property), tolower(water.props("SUPCRT92")))
  if(any(is.na(iprop))) stop(paste("property(s) not available:", paste(property[is.na(iprop)], collapse=" ")))
  # make sure Psat in properties comes with isat=1
  if("psat" %in% tolower(property) & !identical(P, "Psat")) stop("please set P='Psat' to calculate the property Psat")
  # for Psat(T) (1) or T-P (2)
  if(identical(P, "Psat")) iopt <- 1 else iopt <- 2
  if(identical(P, "Psat")) isat <- 1 else isat <- 0
  # input values, gleaned from H2O92D.f and SUP92D.f
  # it, id, ip, ih, itripl, isat, iopt, useLVS, epseqn, icrit
  specs <- c(2, 2, 2, 5, 1, isat, iopt, 1, 4, 0)
  states <- rep(0, 4)
  # initialize the output matrix
  w.out <- matrix(NA, nrow=length(T), ncol=23, byrow=TRUE) 
  err.out <- numeric(length(T))
  rho.out <- numeric(length(T))
  P.out <- numeric(length(T))
  # 20091022 TODO: parallelize this
  Tc <- convert(T, "C")
  for(i in 1:length(T)) {
    states[1] <- Tc[i]
    if(identical(P, "Psat")) states[2] <- 0
    else states[2] <- P[i]
    if(is.na(Tc[i]) | is.na(P[i]) & !identical(P, "Psat")) {
      # if T or P is NA, all properties are NA
      # (NA's are already in w.out)
      P.out[i] <- NA
      rho.out[i] <- NA
    } else {
      # now to the actual calculations
      H2O <- .Fortran("H2O92", as.integer(specs), as.double(states),
        as.double(rep(0, 46)), as.integer(0), PACKAGE="CHNOSZ")
      # errors
      err.out[i] <- H2O[[4]]
      # density of two states
      rho <- H2O[[2]][3]
      rho2 <- H2O[[2]][4]
      if(rho2 > rho) {
        # liquid is denser than vapor
        rho <- rho2 
        inc <- 1  # second state is liquid
      } else inc <- 0  # first state is liquid
      rho.out[i] <- rho
      # 23 properties of the phase in the liquid state
      w <- t(H2O[[3]][seq(1, 45, length.out=23)+inc])
      if(err.out[i]==1) w[1, ] <- NA
      # update the ith row of the output matrix
      w.out[i,] <- w
      # Psat
      if(identical(P, "Psat")) {
        w.P <- H2O[[2]][2]
        w.P[w.P==0] <- NA
        # Psat specifies P=1 below 100 degC
        w.P[w.P < 1] <- 1
        P.out[i] <- w.P
      }
    }
  }
  # convert output to dataframe
  w.out <- as.data.frame(w.out)
  # add names of properties to the output
  names(w.out) <- water.props("SUPCRT92")[1:23]
  # assemble additional properties: V, rho, Psat, E, kT
  if(any(iprop > 23)) {
    mwH2O <- 18.0152 # SUP92.f
    V=mwH2O/rho.out
    rho=rho.out*1000
    Psat=P.out
    E <- V*w.out$alpha
    kT <- V*w.out$beta
    w.out <- cbind(w.out, data.frame(V=V, rho=rho, Psat=Psat, E=E, kT=kT))
  }
  # tell the user about any problems
  if(any(err.out==1)) {
    if(length(T) > 1) plural <- "s" else plural <- ""
    nerr <- length(which(err.out==1))
    if(nerr > 1) plural2 <- "s" else plural2 <- ""
    if(identical(P, "Psat")) msgout(paste("water.SUPCRT92: error", plural2, " calculating ",
      nerr, " of ", length(T), " point", plural, "; for Psat we need 273.16 < T < 647.067 K\n", sep=""))
    else msgout(paste("water.SUPCRT92: error", plural2, " calculating ", nerr,
      " of ", length(T), " point", plural,
      "; T < Tfusion@P, T > 2250 degC, or P > 30kb.\n", sep=""))
      # that last bit is taken from SUP92D.f in SUPCRT92
  }
  # return only the selected properties
  return(w.out[, iprop, drop=FALSE])
}

water.IAPWS95 <- function(property, T=298.15, P=1) {
  # to get the properties of water via IAPWS-95
  msgout(paste("water.IAPWS95: calculating", length(T), "values for"))
  M <- 18.015268 # g mol-1
  v <- function() return(M*1000/my.rho)
  # Psat stuff
  psat <- function() {
    p <- WP02.auxiliary("P.sigma", T)
    p[T < 373.124] <- 0.1
    return(convert(p, "bar"))
  }
  ## thermodynamic properties
  # convert to SUPCRT reference state
  # at the triple point
  # I2S = SUPCRT - IAPWS ( + entropy in G )
  dH <- -68316.76 - 451.75437
  dS <- 16.7123 - 1.581072
  dG <- -56687.71 + 19.64228 - dS * (T - get("thermo")$opt$Tr)
  # does the reference state used for GHS also go here?
  dU <- -67434.5 - 451.3229
  dA <- -55814.06 + 20.07376 - dS * (T - get("thermo")$opt$Tr)
  # calculate pressure from the given T and estimated rho
  p <- function() return(convert(IAPWS95("p", T=T, rho=my.rho), "bar"))
  # convert IAPWS95() (specific, joule) to (molar, cal) 
  s <- function()
    return(convert(IAPWS95('s',T=T,rho=my.rho)$s*M,'cal')+dS) 
  # u (internal energy) is not here because the letter
  # is used to denote one of the Born functions
  # scratch that! let's put u here and call the other one uborn
  u <- function()
    return(convert(IAPWS95('u',T=T,rho=my.rho)$u*M,'cal')+dU)
  a <- function()
    return(convert(IAPWS95('a',T=T,rho=my.rho)$a*M,'cal')+dA)
  h <- function() 
    return(convert(IAPWS95('h',T=T,rho=my.rho)$h*M,'cal')+dH) 
  g <- function() 
    return(convert(IAPWS95('g',T=T,rho=my.rho)$g*M,'cal')+dG) 
  cv <- function() 
    return(convert(IAPWS95('cv',T=T,rho=my.rho)$cv*M,'cal')) 
  cp <- function() 
    return(convert(IAPWS95('cp',T=T,rho=my.rho)$cp*M,'cal')) 
  speed <- function()
    return(IAPWS95('w',T=T,rho=my.rho)$w*100) # to cm/s
  ## electrostatic properties
  diel <- function() return(water.AW90(T=T,rho=my.rho,P=convert(P,'MPa')))
  de.dt <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]
      this.P <- P[i]
      this.rho <- my.rho[i]
      dt <- 0.001; t1 <- this.T-dt; t2 <- this.T+dt
      rho <- rho.IAPWS95(T=c(t1, t2), P=rep(this.P, 2))
      e <- water.AW90(T=c(t1,t2),rho=rho,rep(this.P,2))
      p <- c(p,(e[2]-e[1])/(2*dt))
    }
    return(p)
  }
  de.dp <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]
      this.P <- P[i]
      this.rho <- my.rho[i]
      dp <- 0.001; p1 <- this.P-dp; p2 <- this.P+dp
      rho <- rho.IAPWS95(T=rep(this.T, 2), P=c(p1, p2))
      e <- water.AW90(P=c(p1,p2),rho=rho,T=rep(this.T,2))
      p <- c(p,(e[2]-e[1])/(2*dp))
    }
    return(p)
  }
  ## Born functions
  qborn <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dp <- 0.01; p1 <- this.P-dp; p2 <- this.P+dp
      rho <- rho.IAPWS95(T=rep(this.T, 2), P=c(p1, p2))
      e <- water.AW90(T=rep(this.T,2),rho=rho,P=convert(c(p1,p2),'MPa'))
      #p <- c(p,convert(-(1/e[2]-1/e[1])/(2*dp),'cm3bar'))
      p <- c(p,-(1/e[2]-1/e[1])/(2*dp))
    }
    return(p)
  }
  nborn <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dp <- 0.01; p1 <- this.P-dp; p2 <- this.P+dp
      rho <- rho.IAPWS95(T=rep(this.T, 3), P=c(p1, this.P, p2))
      e <- water.AW90(T=rep(this.T,3),rho=rho,P=convert(c(p1,this.P,p2),'MPa'))
      #p <- c(p,convert(convert((-(1/e[3]-1/e[2])/dp+(1/e[2]-1/e[1])/dp)/dp,'cm3bar'),'cm3bar'))
      p <- c(p,(-(1/e[3]-1/e[2])/dp+(1/e[2]-1/e[1])/dp)/dp)
    }
    return(p)
  }
  yborn <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dt <- 0.001; t1 <- this.T-dt; t2 <- this.T+dt
      rho <- rho.IAPWS95(T=c(t1, t2), P=rep(this.P, 2))
      e <- water.AW90(T=c(t1,t2),rho=rho,P=convert(rep(this.P,2),'MPa'))
      p <- c(p,-(1/e[2]-1/e[1])/(2*dt))
    }
    return(p)
  }
  xborn <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dt <- 0.001; t1 <- this.T-dt; t2 <- this.T+dt
      rho <- rho.IAPWS95(T=c(t1, this.T, t2), P=rep(this.P, 3))
      e <- water.AW90(T=c(t1,this.T,t2),rho=rho,P=convert(rep(this.P,3),'MPa'))
      p <- c(p,(-(1/e[3]-1/e[2])/dt+(1/e[2]-1/e[1])/dt)/dt)
    }
    return(p)
  }
  uborn <- function() {
    p <- numeric()
    for(i in 1:length(T)) {
      this.T <- T[i]; this.P <- P[i]; this.rho <- my.rho[i]
      dt <- 0.001; this.T1 <- this.T - dt; this.T2 <- this.T + dt
      dp <- 0.001; p1 <- this.P-dp; p2 <- this.P+dp
      rho1 <- rho.IAPWS95(T=rep(this.T1, 2), P=c(p1, p2))
      rho2 <- rho.IAPWS95(T=rep(this.T2, 2), P=c(p1, p2))
      e1 <- water.AW90(T=rep(this.T1,2),rho=rho1,P=convert(c(p1,p2),'MPa'))
      e2 <- water.AW90(T=rep(this.T2,2),rho=rho2,P=convert(c(p1,p2),'MPa'))
      #p1 <- convert(-(1/e1[2]-1/e1[1])/(2*dp),'cm3bar')
      #p2 <- convert(-(1/e2[2]-1/e2[1])/(2*dp),'cm3bar')
      p1 <- -(1/e1[2]-1/e1[1])/(2*dp)
      p2 <- -(1/e2[2]-1/e2[1])/(2*dp)
      p <- c(p,(p2-p1)/(2*dt))
    }
    return(p)
  }
  ### main loop; init dataframe output and density holders
  w.out <- NULL
  my.rho <- NULL
  # get densities unless only Psat is requested
  if(!identical(tolower(property), "psat")) {
    # calculate values of P for Psat
    if(identical(P, "Psat")) P <- psat()
    msgout(" rho")
    my.rho <- rho.IAPWS95(T, P, get("thermo")$opt$IAPWS.sat)
    rho <- function() my.rho
  }
  for(i in 1:length(property)) {
    if(tolower(property[i]) %in% c("e", "kt")) {
      # E and kT aren't here yet... set them to NA
      warning("water.IAPWS95: values of ", property[i], " are NA\n", call.=FALSE)
      inew <- rep(NA, length(T))
    } else {
      msgout(paste(" ", property[i], sep=""))
      inew <- get(tolower(property[i]))()
    }
    wnew <- data.frame(inew)
    if(i > 1) w.out <- cbind(w.out, wnew) else w.out <- wnew
  }  
  msgout("\n")
  return(w.out)
}
