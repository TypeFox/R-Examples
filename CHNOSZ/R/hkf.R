# CHNOSZ/hkf.R
# calculate thermodynamic properties using equations of state
# 11/17/03 jmd

hkf <- function(property=NULL,T=298.15,P=1,ghs=NULL,eos=NULL,contrib=c('n','s','o'),
  H2O.PT=NULL,H2O.PrTr=NULL,domega=TRUE) {
  # calculate G, H, S, Cp, V, kT, and/or E using
  # the revised HKF equations of state
  thermo <- get("thermo")
  # constants
  Tr <- thermo$opt$Tr
  Pr <- thermo$opt$Pr
  Theta <- thermo$opt$Theta
  Psi <- thermo$opt$Psi
  # argument handling
  eargs <- eos.args('hkf',property)
  property <- eargs$prop
  props <- eargs$props
  Prop <- eargs$Prop
  domega <- rep(domega,length.out=nrow(ghs))
  # nonsolvation, solvation, and origination contribution
  contribs <- c('n','s','o')
  notcontrib <- ! contrib %in% contribs
  if(TRUE %in% notcontrib)
    stop(paste('argument',c2s(contrib[notcontrib]),'not in',c2s(contribs),'n'))
  # get water properties, if they weren't supplied in arguments (and we want solvation props)
  if('s' %in% contrib) {
    H2O.props <- c("QBorn", "XBorn", "YBorn", "diel")
    # only take these ones if we're in SUPCRT92 compatibility mode
    dosupcrt <- thermo$opt$water != "IAPWS95"
    if(dosupcrt) {
      # (E, daldT, V - for partial derivatives of omega (g function))
      H2O.props <- c(H2O.props,'E','daldT','kT','ZBorn')
    } else {
      # (NBorn, UBorn - for compressibility, expansibility)
      H2O.props <- c(H2O.props,'NBorn','UBorn')
    }
    if(is.null(H2O.PT)) H2O.PT <- water(H2O.props,T=T,P=P)
    if(is.null(H2O.PrTr)) H2O.PrTr <- water(H2O.props,T=thermo$opt$Tr,P=thermo$opt$Pr)
    ZBorn <- -1/H2O.PT$diel
    ZBorn.PrTr <- -1/H2O.PrTr$diel
  }
 # a list to store the result
 x <- list()
 for(k in 1:nrow(ghs)) {
  # loop over each species
  GHS <- ghs[k,]
  EOS <- eos[k,]
  # substitute Cp and V for missing EoS parameters
  # here we assume that the parameters are in the same position as in thermo$obigt
  # put the heat capacity in for c1 if both c1 and c2 are missing
  if(all(is.na(EOS[, 17:18]))) EOS[, 17] <- EOS$Cp
  # put the volume in for a1 if a1, a2, a3 and a4 are missing
  if(all(is.na(EOS[, 13:16]))) EOS[, 13] <- convert(EOS$V, "calories")
  # test for availability of the EoS parameters
  hasEOS <- any(!is.na(EOS[, 13:20]))
  # if at least one of the EoS parameters is available, zero out any NA's in the rest
  if(hasEOS) EOS[, 13:20][, is.na(EOS[, 13:20])] <- 0
  # compute values of omega(P,T) from those of omega(Pr,Tr)
  # using g function etc. (Shock et al., 1992 and others)
  omega <- EOS$omega  # omega.PrTr
  # its derivatives are zero unless the g function kicks in
  dwdP <- dwdT <- d2wdT2 <- numeric(length(T))
  Z <- EOS$Z
  omega.PT <- rep(EOS$omega,length(T))
  if(!is.na(Z)) if(Z != 0) if(domega[k]) if(dosupcrt) {
    # g and f function stuff (Shock et al., 1992; Johnson et al., 1992)
    rhohat <- H2O.PT$rho/1000  # just converting kg/m3 to g/cm3
    g <- gfun(rhohat, convert(T, "C"), P, H2O.PT$alpha, H2O.PT$daldT, H2O.PT$beta)
    # after SUPCRT92/reac92.f
    eta <- 1.66027E5
    reref <- Z^2 / (omega/eta + Z/(3.082 + 0))
    re <- reref + abs(Z) * g$g
    omega.PT <- eta * (Z^2/re - Z/(3.082 + g$g))
    Z3 <- abs(Z^3)/re^2 - Z/(3.082 + g$g)^2
    Z4 <- abs(Z^4)/re^3 - Z/(3.082 + g$g)^3
    dwdP <- (-eta * Z3 * g$dgdP)
    dwdT <- (-eta * Z3 * g$dgdT)
    d2wdT2 <- (2 * eta * Z4 * g$dgdT^2 - eta * Z3 * g$d2gdT2)
  }
  # loop over each property
  w <- NULL
  for(i in 1:length(property)) {
    prop <- property[i]
    # over nonsolvation, solvation, or origination contributions
    hkf.p <- numeric(length(T))
    for(icontrib in contrib) {
      # various contributions to the properties
      if( icontrib=="n") {
        # nonsolvation ghs equations
        if(prop=="h") {
          p.c <- EOS$c1*(T-Tr) - EOS$c2*(1/(T-Theta)-1/(Tr-Theta))
          p.a <- EOS$a1*(P-Pr) + EOS$a2*log((Psi+P)/(Psi+Pr)) + 
            ((2*T-Theta)/(T-Theta)^2)*(EOS$a3*(P-Pr)+EOS$a4*log((Psi+P)/(Psi+Pr)))
          p <- p.c + p.a
        } else if(prop=="s") {
          p.c <- EOS$c1*log(T/Tr) - 
            (EOS$c2/Theta)*( 1/(T-Theta)-1/(Tr-Theta) + 
            log( (Tr*(T-Theta))/(T*(Tr-Theta)) )/Theta )
          p.a <- (T-Theta)^(-2)*(EOS$a3*(P-Pr)+EOS$a4*log((Psi+P)/(Psi+Pr)))
          p <- p.c + p.a
        } else if(prop=="g") {
          p.c <- -EOS$c1*(T*log(T/Tr)-T+Tr) - 
            EOS$c2*( (1/(T-Theta)-1/(Tr-Theta))*((Theta-T)/Theta) - 
            (T/Theta^2)*log((Tr*(T-Theta))/(T*(Tr-Theta))) )
          p.a <- EOS$a1*(P-Pr) + EOS$a2*log((Psi+P)/(Psi+Pr)) + 
            (EOS$a3*(P-Pr) + EOS$a4*log((Psi+P)/(Psi+Pr)))/(T-Theta)
          p <- p.c + p.a
          # at Tr,Pr, if the origination contribution is not NA, ensure the solvation contribution is 0, not NA
          if(!is.na(GHS$G)) p[T==Tr & P==Pr] <- 0
        # nonsolvation cp v kt e equations
        } else if(prop=='cp') {
          p <- EOS$c1 + EOS$c2 * ( T - Theta ) ^ (-2)        
        } else if(prop=='v') {
          p <- convert(EOS$a1,'cm3bar') + 
            convert(EOS$a2,'cm3bar') / ( Psi + P) +
            ( convert(EOS$a3,'cm3bar') + convert(EOS$a4,'cm3bar') / ( Psi + P ) ) / ( T - Theta)
        } else if(prop=='kt') {
          p <- ( convert(EOS$a2,'cm3bar') + 
            convert(EOS$a4,'cm3bar') / (T - Theta) ) * (Psi + P) ^ (-2)
        } else if(prop=='e') {
          p <- convert( - ( EOS$a3 + EOS$a4 / convert((Psi + P),'calories') ) * 
            (T - Theta) ^ (-2),'cm3bar')
        }
      }
      if( icontrib=="s") {
        # solvation ghs equations
        if(prop=="g") {
          p <- -omega.PT*(ZBorn+1) + omega*(ZBorn.PrTr+1) + omega*H2O.PrTr$YBorn*(T-Tr)
          # at Tr,Pr, if the origination contribution is not NA, ensure the solvation contribution is 0, not NA
          if(!is.na(GHS$G)) p[T==Tr & P==Pr] <- 0
        }
        if(prop=="h") 
          p <- -omega.PT*(ZBorn+1) + omega.PT*T*H2O.PT$YBorn + T*(ZBorn+1)*dwdT +
                 omega*(ZBorn.PrTr+1) - omega*Tr*H2O.PrTr$YBorn
        if(prop=="s") 
          p <- omega.PT*H2O.PT$YBorn + (ZBorn+1)*dwdT - omega*H2O.PrTr$YBorn
        # solvation cp v kt e equations
        if(prop=='cp') p <- omega.PT*T*H2O.PT$XBorn + 2*T*H2O.PT$YBorn*dwdT + 
          T*(ZBorn+1)*d2wdT2
        if(prop=='v') p <- -convert(omega.PT,'cm3bar') * 
          H2O.PT$QBorn + convert(dwdP,'cm3bar') * (-ZBorn - 1)
        # WARNING: the partial derivatives of omega are not included here here for kt and e
        # (to do it, see p. 820 of SOJ+92 ... but kt requires d2wdP2 which we don't have yet)
        if(prop=='kt') p <- convert(omega,'cm3bar') * H2O.PT$N
        if(prop=='e') p <- -convert(omega,'cm3bar') * H2O.PT$UBorn
      }
      if( icontrib=='o') {
        # origination ghs equations
        if(prop=='g') {
          p <- GHS$G - GHS$S * (T-Tr)
          # don't inherit NA from GHS$S at Tr
          p[T==Tr] <- GHS$G
        }
        else if(prop=='h') p <- GHS$H
        else if(prop=='s') p <- GHS$S
        # origination eos equations: senseless
        else if(prop %in% tolower(props)) p <- 0 * T
      }
      # accumulate the contribution
      hkf.p <- hkf.p + p
    }
    wnew <- data.frame(hkf.p)
    if(i>1) w <- cbind(w,wnew) else w <- wnew
  }
  colnames(w) <- Prop
  x[[k]] <- w
 }
 return(x)
}

gfun <- function(rhohat, Tc, P, alpha, daldT, beta) {
  ## g and f functions for describing effective electrostatic radii of ions
  ## split from hkf() 20120123 jmd      
  ## based on equations in
  ## Shock EL, Oelkers EH, Johnson JW, Sverjensky DA, Helgeson HC, 1992
  ## Calculation of the Thermodynamic Properties of Aqueous Species at High Pressures 
  ## and Temperatures: Effective Electrostatic Radii, Dissociation Constants and 
  ## Standard Partial Molal Properties to 1000 degrees C and 5 kbar
  ## J. Chem. Soc. Faraday Trans., 88(6), 803-826  doi:10.1039/FT9928800803
  # rhohat - density of water in g/cm3
  # Tc - temperature in degrees Celsius
  # P - pressure in bars
  # start with an output list of zeros
  out0 <- numeric(length(rhohat))
  out <- list(g=out0, dgdT=out0, d2gdT2=out0, dgdP=out0)
  # only rhohat less than 1 will give results other than zero
  idoit <- rhohat < 1 & !is.na(rhohat)
  rhohat <- rhohat[idoit]
  Tc <- Tc[idoit]
  P <- P[idoit]
  alpha <- alpha[idoit]
  daldT <- daldT[idoit]
  beta <- beta[idoit]
  # eta in Eq. 1
  eta <- 1.66027E5
  # Table 3
  ag1 <- -2.037662
  ag2 <- 5.747000E-3
  ag3 <- -6.557892E-6
  bg1 <- 6.107361
  bg2 <- -1.074377E-2
  bg3 <- 1.268348E-5
  # Eq. 25
  ag <- ag1 + ag2 * Tc + ag3 * Tc ^ 2
  # Eq. 26
  bg <- bg1 + bg2 * Tc + bg3 * Tc ^ 2
  # Eq. 24
  g <- ag * (1 - rhohat) ^ bg
  # Table 4
  af1 <- 0.3666666E2
  af2 <- -0.1504956E-9
  af3 <- 0.5017997E-13
  # Eq. 33
  f <- 
    ( ((Tc - 155) / 300) ^ 4.8 + af1 * ((Tc - 155) / 300) ^ 16 ) *
    ( af2 * (1000 - P) ^ 3 + af3 * (1000 - P) ^ 4 ) 
  # limits of the f function (region II of Fig. 6)
  ifg <- Tc > 155 & P < 1000 & Tc < 355
  # in case any T or P are NA
  ifg <- ifg & !is.na(ifg)
  # Eq. 32
  g[ifg] <- g[ifg] - f[ifg]
  ## now we have g at P, T
  ## the rest is to get its partial derivatives with pressure and temperature
  ## after Johnson et al., 1992
  # alpha - coefficient of isobaric expansivity (K^-1)
  # daldT - temperature derivative of coefficient of isobaric expansivity (K^-2)
  # beta - coefficient of isothermal compressibility (bar^-1)
  # Eqn. 76
  d2fdT2 <- (0.0608/300*((Tc-155)/300)^2.8 + af1/375*((Tc-155)/300)^14) * (af2*(1000-P)^3 + af3*(1000-P)^4)
  # Eqn. 75
  dfdT <- (0.016*((Tc-155)/300)^3.8 + 16*af1/300*((Tc-155)/300)^15) * (af2*(1000-P)^3 + af3*(1000-P)^4)
  # Eqn. 74
  dfdP <- -(((Tc-155)/300)^4.8 + af1*((Tc-155)/300)^16) * (3*af2*(1000-P)^2 + 4*af3*(1000-P)^3)
  d2bdT2 <- 2 * bg3  # Eqn. 73
  d2adT2 <- 2 * ag3  # Eqn. 72
  dbdT <- bg2 + 2*bg3*Tc  # Eqn. 71
  dadT <- ag2 + 2*ag3*Tc  # Eqn. 70
  # Eqn. 69
  dgadT <- bg*rhohat*alpha*(1-rhohat)^(bg-1) + log(1-rhohat)*g/ag*dbdT  
  D <- rhohat
  # transcribed from SUPCRT92/reac92.f
  dDdT <- -D * alpha
  dDdP <- D * beta
  dDdTT <- -D * (daldT - alpha^2)
  Db <- (1-D)^bg
  dDbdT <- -bg*(1-D)^(bg-1)*dDdT + log(1-D)*Db*dbdT
  dDbdTT <- -(bg*(1-D)^(bg-1)*dDdTT + (1-D)^(bg-1)*dDdT*dbdT + 
    bg*dDdT*(-(bg-1)*(1-D)^(bg-2)*dDdT + log(1-D)*(1-D)^(bg-1)*dbdT)) +
    log(1-D)*(1-D)^bg*d2bdT2 - (1-D)^bg*dbdT*dDdT/(1-D) + log(1-D)*dbdT*dDbdT
  d2gdT2 <- ag*dDbdTT + 2*dDbdT*dadT + Db*d2adT2
  d2gdT2[ifg] <- d2gdT2[ifg] - d2fdT2[ifg]
  dgdT <- g/ag*dadT + ag*dgadT  # Eqn. 67
  dgdT[ifg] <- dgdT[ifg] - dfdT[ifg]
  dgdP <- -bg*rhohat*beta*g*(1-rhohat)^-1  # Eqn. 66
  dgdP[ifg] <- dgdP[ifg] - dfdP[ifg]
  # phew! done with those derivatives
  # put the results in their right place (where rhohat < 1)
  out$g[idoit] <- g
  out$dgdT[idoit] <- dgdT
  out$d2gdT2[idoit] <- d2gdT2
  out$dgdP[idoit] <- dgdP
  return(out)
}

