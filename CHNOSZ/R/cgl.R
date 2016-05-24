# CHNOSZ/cgl.R
# calculate standard thermodynamic properties of non-aqueous species
# 20060729 jmd

cgl <- function(property=NULL,T=298.15,P=1,ghs=NULL,eos=NULL) {
  # calculate properties of crystalline, liquid (except H2O) and gas species
  # argument handling
  thermo <- get("thermo")
  Tr <- thermo$opt$Tr
  Pr <- thermo$opt$Pr
  eargs <- eos.args('mk',property=property)
  prop <- eargs$prop
  props <- eargs$props
  Prop <- eargs$Prop

  # a list - the result
  x <- list()
  for(k in 1:nrow(ghs)) {
    # loop over each species
    GHS <- ghs[k,]
    EOS <- eos[k,]
    w <- NULL
    for(i in 1:length(prop)) {
      property <- prop[i]
      # a test for availability of the EoS parameters
      # here we assume that the parameters are in the same position as in thermo$obigt
      # leave T transition (in 20th column) alone
      hasEOS <- any(!is.na(EOS[, 13:19]))
      # if at least one of the EoS parameters is available, zero out any NA's in the rest
      if(hasEOS) EOS[, 13:19][, is.na(EOS[, 13:19])] <- 0
      # equations for lambda adapted from HOK+98
      if(property=='cp') {
        # use constant Cp if the EoS parameters are not available
        if(!hasEOS) p <- EOS$Cp
        else p <- EOS$a + EOS$b*T + EOS$c*T^-2 + EOS$d*T^-0.5 + EOS$e*T^2 + EOS$f*T^EOS$lambda
      } else if(property=='v') {
        p <- rep(EOS$V,length(T))
      } else if(property %in% c('e','kt')) {
        p <- rep(NA,length(T))
        warning('cgl: E and/or kT of cr, gas and/or liq species are NA.')
      } else if(property=='g') {
        # use constant Cp if the EoS parameters are not available
        if(!hasEOS) p <- GHS$G + EOS$Cp*(T-Tr-T*log(T/Tr)) else
        # Gibbs energy integral: the value at Tref plus heat capacity terms
        p <-   GHS$G + EOS$a*(T-Tr-T*log(T/Tr)) - 
               EOS$b*(T-Tr)^2/2 - EOS$c*(1/T + T/Tr^2 - 2/Tr)/2 -
               EOS$d*(T^0.5-0.5*T*Tr^-0.5-0.5*Tr^0.5)/-0.25 -
               EOS$e*(T^3-3*T*Tr^2+2*Tr^3)/6
        # entropy and volume terms
        if(!is.na(GHS$S)) p <- p - GHS$S*(T-Tr)
        if(!is.na(EOS$V)) p <- p + convert(EOS$V*(P-Pr),'calories')
        # use additional heat capacity term if it's defined
        if(!is.na(EOS$f) & !is.na(EOS$lambda)) if(EOS$f!=0) {
          if(EOS$lambda== -1) p <- p + EOS$f*(log(T/Tr)-T*(1/Tr-1/T))
          else p <- p + EOS$f * ( T^(EOS$lambda+1) - (EOS$lambda+1)*T*Tr^EOS$lambda + 
            EOS$lambda*Tr^(EOS$lambda+1) ) / ( EOS$lambda*(EOS$lambda+1) ) 
        }
      } else if(property=='h') { 
        # use constant Cp if the EoS parameters are not available
        if(!hasEOS) p <- EOS$Cp*(T-Tr) else
        p <- EOS$a*(T-Tr) + EOS$b*(T^2-Tr^2)/2 +
             EOS$c*(1/T-1/Tr)/-1 + EOS$d*(T^0.5-Tr^0.5)/0.5 + 
             EOS$e*(T^3-Tr^3)/3 
             # SUPCRT seems to ignore this term? ... 20070802
             # + convert(EOS$V*(P-Pr),'calories')
        p <- GHS$H + p
        if(!is.na(EOS$f) & !is.na(EOS$lambda)) if(EOS$f!=0) {
           if(EOS$lambda == -1) p <- p + EOS$f*log(T/Tr) 
           else p <- p - EOS$f * ( T^(EOS$lambda+1) - Tr^(EOS$lambda+1) ) / (EOS$lambda+1)
        }
      } else if(property=='s') {
        # use constant Cp if the EoS parameters are not available
        if(!hasEOS) p <- EOS$Cp*log(T/Tr) else
        p <- EOS$a*log(T/Tr) + EOS$b*(T-Tr) + 
             EOS$c*(T^(-2)-Tr^(-2))/(-2) + EOS$e*(T^2-Tr^2)/2 + 
             EOS$d*(T^-0.5-Tr^-0.5)/-0.5
        p <- GHS$S + p
        if(!is.na(EOS$f) & !is.na(EOS$lambda)) if(EOS$f!=0) {
          p <- p + EOS$f*(T^EOS$lambda-Tr^EOS$lambda)/EOS$lambda
        }
      }
      wnew <- data.frame(p)
      if(i>1) w <- cbind(w,wnew) else w <- wnew
    }
  colnames(w) <- Prop
  x[[k]] <- w
 }
 return(x)
}

