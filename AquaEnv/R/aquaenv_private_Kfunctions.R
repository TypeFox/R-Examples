#######################################################################################
# local functions
######################################################################################
lnK <- function(A, B, C, D, E, T)        
  {
    lnK <- A + (B/T) + C*log(T) + D*T + E*T^2
  }

logK <- function(A, B, C, D, E, T)        
  {
    logK <- A + (B/T) + C*log10(T) + D*T + E*T^2
  }

Sterms <- function(S)
  {
    sqrtS  <- sqrt(S)
    Sterms <- list(S^2, sqrtS, S*sqrtS)
  }

Iterms <- function(S)
  {
    I      <- I(S)
    sqrtI  <- sqrt(I)
    Iterms <- list(I, I^2, sqrtI, I*sqrtI)
  }

deltaPlnK <- function(T, p, coeff)
  {
    t         <- T + PhysChemConst$absZero 
    t2        <- t^2
    deltaV    <- coeff[[1]] + coeff[[2]]*t + coeff[[3]]*t2
    deltaK    <- (coeff[[4]] + coeff[[5]]*t + coeff[[6]]*t2)/1000 #the 1000 is correction from Lewis1998
    deltaPlnK <- -(deltaV/(PhysChemConst$R*T))*p + (0.5*(deltaK/(PhysChemConst$R*T)))*(p^2)  
  }

splitS_K_CO2 <- function(T)
  {
    res <- c()
    for (te in T)
      {
        lnK1 <- function(T,S)
          {
            Sterms <- Sterms(S)
            
            A <- 2.83655    - 0.20760841*Sterms[[2]] + 0.08468345*S - 0.00654208*Sterms[[3]]
            B <- -2307.1266 -     4.0484*Sterms[[2]]
            C <- -1.5529413
            
            lnK1 <- lnK(A, B, C, 0, 0, T)    
          }
        
        lnK2 <- function(T,S)
          {
            Sterms <- Sterms(S)
            
            A <- 290.9097  -  228.39774*Sterms[[2]] +   54.20871*S -  3.969101*Sterms[[3]] - 0.00258768*Sterms[[1]]
            B <- -14554.21 + 9714.36839*Sterms[[2]] - 2310.48919*S + 170.22169*Sterms[[3]]
            C <- -45.0575  +  34.485796*Sterms[[2]] -    8.19515*S +   0.60367*Sterms[[3]]
            
            lnK2 <- lnK(A, B, C, 0, 0, T)    
          }
        res <- c(res, uniroot(function(x) lnK1(te,x)-lnK2(te,x), c(3,7))$root)
      }
    return(res)
  }

splitS_K_HCO3 <- function(T)
  {
    res <- c()
    for (te in T)
      {
        lnK1 <- function(T,S)
          {
            Sterms <- Sterms(S)
            
            A <- -9.226508  - 0.106901773*Sterms[[2]] + 0.1130822*S - 0.00846934*Sterms[[3]] 
            B <- -3351.6106 -     23.9722*Sterms[[2]] 
            C <- -0.2005743 
            
            lnK1 <- lnK(A, B, C, 0, 0, T)    
          }
        
        lnK2 <- function(T,S)
          {
            Sterms <- Sterms(S)
            
            A <- 207.6548  -  167.69908*Sterms[[2]] +   39.75854*S -   2.892532*Sterms[[3]] - 0.00613142*Sterms[[1]]
            B <- -11843.79 + 6551.35253*Sterms[[2]] - 1566.13883*S + 116.270079*Sterms[[3]]
            C <- -33.6485  +  25.928788*Sterms[[2]] -   6.171951*S + 0.45788501*Sterms[[3]]
            
            lnK2 <- lnK(A, B, C, 0, 0, T)    
          }
        res <- c(res, uniroot(function(x) lnK1(te,x)-lnK2(te,x), c(3,7))$root)
      }
    return(res)
  }

att <- function(K)
  {
    attr(K, "unit")     <- "mol/kg-soln"
    attr(K, "pH scale") <- "free"
    att                 <- K
  }
######################################################################################





######################################################################################
# Henry's constants
######################################################################################
K0_CO2 <- function(S, t)           
  {
    T <- T(t)

    A <- 0.023517*S - 167.81077
    B <- 9345.17
    C <- 23.3585
    D <- -2.3656e-4*S
    E <- 4.7036e-7*S
    
    K0_CO2 <- exp(lnK(A, B, C, D, E, T))
    attr(K0_CO2, "unit") <- "mol/(kg-soln*atm)"
    
    return (K0_CO2)
  }

K0_O2 <- function(S, t)            
  {
    T <- T(t) 
    
    A <- -846.9978 - 0.037362*S
    B <- 25559.07 
    C <- 146.4813
    D <- -0.22204 + 0.00016504*S
    E <- -2.0564e-7*S
    
    K0_O2 <- exp(lnK(A, B, C, D, E, T))
    K0_O2 <- K0_O2 * PhysChemConst$uMolToMol
    attr(K0_O2, "unit") <- "mol/(kg-soln*atm)"

    return(K0_O2)
  }
######################################################################################




######################################################################################
# ion product of water
######################################################################################
K_W <- function(S, t, p=0, SumH2SO4=NULL, SumHF=NULL, khf="dickson", khso4="dickson")
  {
    T <- T(t) 
    Sterms <- Sterms(S)
    
    A <- 148.9652  -  5.977*Sterms[[2]] - 0.01615*S
    B <- -13847.26 + 118.67*Sterms[[2]]
    C <- -23.6521  + 1.0495*Sterms[[2]]
    D <- 0
    E <- 0
    
    K_W <- exp(lnK(A, B, C, D, E, T) +
               log(scaleconvert(S, t, p=0, SumH2SO4, SumHF, khf, khso4)$tot2sws) + 
               deltaPlnK(T, p, DeltaPcoeffs$K_W) +
               log(scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$sws2free))
    attr(K_W, "unit")     <- "(mol/kg-soln)^2"
    attr(K_W, "pH scale") <- "free" 
    
    return(K_W)
  }
######################################################################################




######################################################################################
# acid dissociation constants
######################################################################################
K_HSO4 <- function(S, t, p=0, khso4="dickson")           
  {
    if (khso4 == "khoo")
      {
        return(K_HSO4_khoo(S=S, t=t, p=p))
      }
    else
      {
        T <- T(t)
        Iterms <- Iterms(S)
        
        A <-  324.57*Iterms[[3]] - 771.54*Iterms[[1]] + 141.328  
        B <-   35474*Iterms[[1]] +   1776*Iterms[[2]] - 13856*Iterms[[3]] - 2698*Iterms[[4]] - 4276.1
        C <- 114.723*Iterms[[1]] - 47.986*Iterms[[3]] - 23.093
        D <- 0
        E <- 0
        
        K_HSO4 <- exp(lnK(A, B, C, D, E, T) + deltaPlnK(T, p, DeltaPcoeffs$K_HSO4) + log(molal2molin(S)))
        
        return (eval(att(K_HSO4)))
      }
  }

K_HSO4_khoo <- function(S, t, p=0)           
  {
    T <- T(t)
        
    A <- 6.3451 + 0.5208*sqrt(I(S))
    B <- -647.59 
    C <- 0
    D <- -0.019085
    E <- 0
    
    K_HSO4_khoo <- exp(lnK(A, B, C, D, E, T) + deltaPlnK(T, p, DeltaPcoeffs$K_HSO4) + log(molal2molin(S)))

    return (eval(att(K_HSO4_khoo)))
  }

K_HF <- function(S, t, p=0, SumH2SO4=NULL, SumHF=NULL, khf="dickson", khso4="dickson")
  {
    if (khf == "perez")
      {
        return(K_HF_perez(S=S, t=t, p=p, SumH2SO4=SumH2SO4, SumHF=SumHF, khso4=khso4))
      }
    else
      {
        T <- T(t)
        
        A <- 1.525 * sqrt(I(S)) - 12.641
        B <- 1590.2
        C <- 0
        D <- 0
        E <- 0
        
        K_HF <- exp(lnK(A, B, C, D, E, T) + deltaPlnK(T, p, DeltaPcoeffs$K_HF) + log(molal2molin(S)))
        
        return(eval(att(K_HF)))
      }
  }

K_HF_perez <- function(S, t, p=0, SumH2SO4=NULL, SumHF=NULL, khso4="dickson")
  {
    T <- T(t)
    
    A <- -9.68  + 0.111 * sqrt(S) 
    B <- 874
    C <- 0
    D <- 0
    E <- 0

    # To convert from the total to the free scale, only K_HSO4 is needed. That's why we do not need to call "scaleconvert" with
    # the option khf="perez" here. Otherwise the cat would also bite its own tail
    K_HF_perez <- exp(lnK(A, B, C, D, E, T) + deltaPlnK(T, p, DeltaPcoeffs$K_HF) +
                      log(scaleconvert(S, t, p=0, SumH2SO4, SumHF, khso4=khso4)$tot2free))        # CODE of CO2Sys matlab version:
                                                                                     # one can assume that the pressure corrections for K_HF and K_HSO4 given in Millero 1995 are
                                                                                     # valid for constants on the FREE pH scale only, that's why we convert the constant
                                                                                     # first from the total to the free scale WITHOUT pressure corrected conversion factors
                                                                                     # and then pressure correct (summands are switched here). If one would now want to
                                                                                     # convert this obtained constant to another scale, it would need to be done WITH
                                                                                     # pressure corrected conversion factors (as implemented in "convert")
    return(eval(att(K_HF_perez)))
  }

K_CO2 <- function(S, t, p=0, SumH2SO4=NULL, SumHF=NULL, k1k2="roy", khf="dickson", khso4="dickson") # Roy1993
  {
    if (k1k2 == "lueker")
      {
        return(K_CO2_lueker(S=S, t=t, p=p, SumH2SO4=SumH2SO4, SumHF=SumHF, khf=khf, khso4=khso4))
      }
    else if (k1k2 == "millero")
      {
        return(K_CO2_millero(S=S, t=t, p=p, SumH2SO4=SumH2SO4, SumHF=SumHF, khf=khf, khso4=khso4))
      }
    else
      { 
        T     <- T(t)
        splitS <- splitS_K_CO2(T)
        K_CO2  <- c()
        
        for (s in 1:length(S))
          {
            Sterms <- Sterms(S[[s]])
            
            for (te in 1:length(t))
              {
                if (S[[s]] > splitS[[te]])
                  {
                    A <- 2.83655    - 0.20760841*Sterms[[2]] + 0.08468345*S[[s]] - 0.00654208*Sterms[[3]]
                    B <- -2307.1266 -     4.0484*Sterms[[2]]
                    C <- -1.5529413
                    D <- 0
                    E <- 0
                  }
                else
                  {
                    A <- 290.9097  -  228.39774*Sterms[[2]] +   54.20871*S[[s]] -  3.969101*Sterms[[3]] - 0.00258768*Sterms[[1]]
                    B <- -14554.21 + 9714.36839*Sterms[[2]] - 2310.48919*S[[s]] + 170.22169*Sterms[[3]]
                    C <- -45.0575  +  34.485796*Sterms[[2]] -    8.19515*S[[s]] +   0.60367*Sterms[[3]]
                    D <- 0
                    E <- 0
                  }
              }
            
            K_CO2 <- c(K_CO2, exp(lnK(A, B, C, D, E, T) +
                                  log(scaleconvert(S[[s]], t, p=0, SumH2SO4[[s]], SumHF[[s]], khf=khf, khso4=khso4)$tot2sws) +    #Lewis2008: first convert to SWS with not pressure corrected conversion
                                  deltaPlnK(T, p, DeltaPcoeffs$K_CO2) +                                                           #           then pressure correct (inferred from Millero1995: valid on SWS)
                                  log(scaleconvert(S[[s]], t, p, SumH2SO4[[s]], SumHF[[s]], khf=khf, khso4=khso4)$sws2free) +     #           then convert to destiny scale with pressure corrected conversion
                                  log(molal2molin(S[[s]]))))
          }
        return(eval(att(K_CO2)))
      }
  }

K_CO2_lueker <- function(S, t, p=0, SumH2SO4=NULL, SumHF=NULL, khf="dickson", khso4="dickson") # Lueker2000
  {
    T      <- T(t)
    Sterms <- Sterms(S)

    A <- 61.2172 + 0.011555*S - 0.0001152*Sterms[[1]]
    B <- -3633.86
    C <- -9.67770
    D <- 0
    E <- 0
            
    K_CO2_lueker <- 10^(lnK(A, B, C, D, E, T)) * exp(log(scaleconvert(S, t, p=0, SumH2SO4, SumHF, khf=khf, khso4=khso4)$tot2sws) +
                                                     deltaPlnK(T, p, DeltaPcoeffs$K_CO2) +
                                                     log(scaleconvert(S, t, p, SumH2SO4, SumHF, khf=khf, khso4=khso4)$sws2free))
    
    return(eval(att(K_CO2_lueker)))
  }

K_CO2_millero <- function(S, t, p=0, SumH2SO4=NULL, SumHF=NULL, khf="dickson", khso4="dickson") # Millero2006
  {
    T      <- T(t)
    Sterms <- Sterms(S) 

    A <- 126.34048 - 0.0331*S + 0.0000533*Sterms[[1]] - 13.4191*Sterms[[2]]
    B <- -6320.813 +  6.103*S                         + 530.123*Sterms[[2]]
    C <- -19.568224                                   + 2.06950*Sterms[[2]]
    D <- 0
    E <- 0
            
    K_CO2_millero <- 10^(lnK(A, B, C, D, E, T)) * exp(deltaPlnK(T, p, DeltaPcoeffs$K_CO2) + log(scaleconvert(S, t, p, SumH2SO4, SumHF, khf=khf, khso4=khso4)$sws2free))
    
    return(eval(att(K_CO2_millero)))
  }

K_HCO3 <- function(S, t, p=0, SumH2SO4=NULL, SumHF=NULL, k1k2="roy", khf="dickson", khso4="dickson") # Roy1993
  {
    if (k1k2 == "lueker")
      {
        return(K_HCO3_lueker(S=S, t=t, p=p, SumH2SO4=SumH2SO4, SumHF=SumHF, khf=khf, khso4=khso4))
      }
    else if (k1k2 == "millero")
      {
        return(K_HCO3_millero(S=S, t=t, p=p, SumH2SO4=SumH2SO4, SumHF=SumHF, khf=khf, khso4=khso4))
      }
    else
      {
        T     <- T(t)
        splitS <- splitS_K_HCO3(T)
        K_HCO3 <- c()
        
        for (s in 1:length(S))
          {
            Sterms <- Sterms(S[[s]])
            
            for (te in 1:length(t))
              {
                if (S[[s]] > splitS[[te]])
                  {
                    A <- -9.226508  - 0.106901773*Sterms[[2]] + 0.1130822*S[[s]] - 0.00846934*Sterms[[3]] 
                    B <- -3351.6106 -     23.9722*Sterms[[2]] 
                    C <- -0.2005743 
                    D <- 0
                    E <- 0
                  }
                else
                  {
                    A <- 207.6548  -  167.69908*Sterms[[2]] +   39.75854*S[[s]] -   2.892532*Sterms[[3]] - 0.00613142*Sterms[[1]]
                    B <- -11843.79 + 6551.35253*Sterms[[2]] - 1566.13883*S[[s]] + 116.270079*Sterms[[3]]
                    C <- -33.6485  +  25.928788*Sterms[[2]] -   6.171951*S[[s]] + 0.45788501*Sterms[[3]]
                    D <- 0
                    E <- 0
                  }
              }
            
            K_HCO3 <- c(K_HCO3, exp(lnK(A, B, C, D, E, T) +
                                    log(scaleconvert(S[[s]], t, p=0, SumH2SO4[[s]], SumHF[[s]], khf=khf, khso4=khso4)$tot2sws) +
                                    deltaPlnK(T, p, DeltaPcoeffs$K_HCO3) +
                                    log(scaleconvert(S[[s]], t, p, SumH2SO4[[s]], SumHF[[s]], khf=khf, khso4=khso4)$sws2free) +
                                    log(molal2molin(S[[s]]))))
          }
        
        return(eval(att(K_HCO3)))
      }
  }

K_HCO3_lueker <- function(S, t, p=0, SumH2SO4=NULL, SumHF=NULL, khf="dickson", khso4="dickson") # Lueker2000
  {
    T      <- T(t)
    Sterms <- Sterms(S)

    A <- -25.9290 + 0.01781*S -0.0001122*Sterms[[1]]
    B <- -471.78
    C <- 3.16967
    D <- 0
    E <- 0
            
    K_HCO3_lueker <- 10^(lnK(A, B, C, D, E, T)) * exp(log(scaleconvert(S, t, p=0, SumH2SO4, SumHF, khf=khf, khso4=khso4)$tot2sws) +
                                                      deltaPlnK(T, p, DeltaPcoeffs$K_HCO3) +
                                                      log(scaleconvert(S, t, p, SumH2SO4, SumHF, khf=khf, khso4=khso4)$sws2free))
    
    return(eval(att(K_HCO3_lueker)))
  }

K_HCO3_millero <- function(S, t, p=0, SumH2SO4=NULL, SumHF=NULL, khf="dickson", khso4="dickson") # Millero2006
  {
    T      <- T(t)
    Sterms <- Sterms(S) 

    A <- 90.18333  - 0.1248*S + 0.0003687*Sterms[[1]] - 21.0894*Sterms[[2]]
    B <- -5143.692 + 20.051*S                         + 772.483*Sterms[[2]]
    C <- -14.613358                                   +  3.3336*Sterms[[2]]
    D <- 0
    E <- 0
            
    K_HCO3_millero <- 10^(lnK(A, B, C, D, E, T)) * exp(deltaPlnK(T, p, DeltaPcoeffs$K_HCO3) + log(scaleconvert(S, t, p, SumH2SO4, SumHF, khf=khf, khso4=khso4)$sws2free))
    
    return(eval(att(K_HCO3_millero)))
  }

K_BOH3 <- function(S, t, p=0, SumH2SO4=NULL, SumHF=NULL, khf="dickson", khso4="dickson")
  {
    T <- T(t)
    Sterms <- Sterms(S)
    
    A <- 148.0248 + 137.1942*Sterms[[2]] + 1.62142*S
    B <- -8966.90 -  2890.53*Sterms[[2]] -  77.942*S + 1.728*Sterms[[3]] - 0.0996*Sterms[[1]]
    C <- -24.4344 -   25.085*Sterms[[2]] -  0.2474*S
    D <-            0.053105*Sterms[[2]]
    E <- 0
    
    K_BOH3 <- exp(lnK(A, B, C, D, E, T) +
                  log(scaleconvert(S, t, p=0, SumH2SO4, SumHF, khf=khf, khso4=khso4)$tot2sws) +
                  deltaPlnK(T, p, DeltaPcoeffs$K_BOH3) +
                  log(scaleconvert(S, t, p, SumH2SO4, SumHF, khf=khf, khso4=khso4)$sws2free))
    
    return(eval(att(K_BOH3)))
  }

K_NH4 <- function(S, t, p=0, SumH2SO4=NULL, SumHF=NULL, khf="dickson", khso4="dickson")
  {
    T <- T(t)
    Sterms <- Sterms(S)
    
    A <- -0.25444 +  0.46532*Sterms[[2]] - 0.01992*S
    B <- -6285.33 - 123.7184*Sterms[[2]] + 3.17556*S
    C <- 0
    D <- 0.0001635
    E <- 0
    
    K_NH4 <- exp(lnK(A, B, C, D, E, T) + deltaPlnK(T, p, DeltaPcoeffs$K_NH4) + log(scaleconvert(S, t, p, SumH2SO4, SumHF, khf=khf, khso4=khso4)$sws2free))
    
    return(eval(att(K_NH4)))
  }

K_H2S <- function(S, t, p=0, SumH2SO4=NULL, SumHF=NULL, khf="dickson", khso4="dickson")
  {
    T <- T(t)
    Sterms <- Sterms(S)
    
    A <- 225.838  + 0.3449*Sterms[[2]] - 0.0274*S
    B <- -13275.3
    C <- -34.6435
    D <- 0
    E <- 0
    
    K_H2S <- exp(lnK(A, B, C, D, E, T) +
                 log(scaleconvert(S, t, p=0, SumH2SO4, SumHF, khf=khf, khso4=khso4)$tot2sws) +
                 deltaPlnK(T, p, DeltaPcoeffs$K_H2S) +
                 log(scaleconvert(S, t, p, SumH2SO4, SumHF, khf=khf, khso4=khso4)$sws2free))
        
    return(eval(att(K_H2S)))
  }

K_H3PO4 <- function(S, t, p=0, SumH2SO4=NULL, SumHF=NULL, khf="dickson", khso4="dickson")
  {
    T <- T(t)
    Sterms <- Sterms(S)
    
    A <- 115.525   + 0.69171*Sterms[[2]] - 0.01844*S
    B <- -4576.752 - 106.736*Sterms[[2]] - 0.65643*S
    C <- -18.453
    D <- 0
    E <- 0
    
    K_H3PO4 <- exp(lnK(A, B, C, D, E, T) +
                   log(scaleconvert(S, t, p=0, SumH2SO4, SumHF, khf=khf, khso4=khso4)$tot2sws) +
                   deltaPlnK(T, p, DeltaPcoeffs$K_H3PO4) +
                   log(scaleconvert(S, t, p, SumH2SO4, SumHF, khf=khf, khso4=khso4)$sws2free))
    
    return(eval(att(K_H3PO4)))
  }

K_H2PO4 <- function(S, t, p=0, SumH2SO4=NULL, SumHF=NULL, khf="dickson", khso4="dickson") 
  {
    T <- T(t)
    Sterms <- Sterms(S)
    
    A <- 172.0883  +  1.3566*Sterms[[2]] - 0.05778*S
    B <- -8814.715 - 160.340*Sterms[[2]] + 0.37335*S
    C <- -27.927
    D <- 0
    E <- 0
    
    K_H2PO4 <- exp(lnK(A, B, C, D, E, T) +
                   log(scaleconvert(S, t, p=0, SumH2SO4, SumHF, khf=khf, khso4=khso4)$tot2sws) +
                   deltaPlnK(T, p, DeltaPcoeffs$K_H2PO4) +
                   log(scaleconvert(S, t, p, SumH2SO4, SumHF, khf=khf, khso4=khso4)$sws2free))
    
    return(eval(att(K_H2PO4)))
  }

K_HPO4 <- function(S, t, p=0, SumH2SO4=NULL, SumHF=NULL, khf="dickson", khso4="dickson")
  {
    T <- T(t)
    Sterms <- Sterms(S)
    
    A <- -18.141  +  2.81197*Sterms[[2]] -  0.09984*S
    B <- -3070.75 + 17.27039*Sterms[[2]] - 44.99486*S
    C <- 0
    D <- 0
    E <- 0 
    
    K_HPO4 <- exp(lnK(A, B, C, D, E, T) +
                  log(scaleconvert(S, t, p=0, SumH2SO4, SumHF, khf=khf, khso4=khso4)$tot2sws) +
                  deltaPlnK(T, p, DeltaPcoeffs$K_HPO4) +
                  log(scaleconvert(S, t, p, SumH2SO4, SumHF, khf=khf, khso4=khso4)$sws2free))
    
    return(eval(att(K_HPO4)))
  }

K_SiOH4 <- function(S, t, p=0, SumH2SO4=NULL, SumHF=NULL, khf="dickson", khso4="dickson")
  {
    T <- T(t)
    Iterms <- Iterms(S)
    
    A <- 117.385 + 3.5913*Iterms[[3]] - 1.5998*Iterms[[1]] + 0.07871*Iterms[[2]]
    B <- -8904.2 - 458.79*Iterms[[3]] + 188.74*Iterms[[1]] - 12.1652*Iterms[[2]]
    C <- -19.334
    D <- 0
    E <- 0
    
    K_SiOH4 <- exp(lnK(A, B, C, D, E, T) +
                   log(scaleconvert(S, t, p=0, SumH2SO4, SumHF, khf=khf, khso4=khso4)$tot2sws) + 
                   deltaPlnK(T, p, DeltaPcoeffs$K_SiOH4) +
                   log(scaleconvert(S, t, p, SumH2SO4, SumHF, khf=khf, khso4=khso4)$sws2free) +
                   log(molal2molin(S)))
    
    return(eval(att(K_SiOH4)))
  }

K_SiOOH3 <- function(S, t, p=0, SumH2SO4=NULL, SumHF=NULL, khf="dickson", khso4="dickson")
  {
    T <- T(t)

    A <- 8.96
    B <- -4465.18
    C <- 0
    D <- 0.021952
    E <- 0
    
    K_SiOOH3 <- exp(lnK(A, B, C, D, E, T) +
                    log(scaleconvert(S, t, p=0, SumH2SO4, SumHF, khf=khf, khso4=khso4)$tot2sws) + 
                    deltaPlnK(T, p, DeltaPcoeffs$K_SiOOH3) +
                    log(scaleconvert(S, t, p, SumH2SO4, SumHF, khf=khf, khso4=khso4)$sws2free) +
                    log(molal2molin(S)))
    
    return(eval(att(K_SiOOH3)))
  }
######################################################################################



######################################################################################
# solubility products
######################################################################################
Ksp_calcite <- function(S, t, p=0) 
  {
    T <- T(t)
    Sterms <- Sterms(S)
    
    A <- -171.9065 -   0.77712*Sterms[[2]] - 0.07711*S + 0.0041249*Sterms[[3]]
    B <- 2839.319  +    178.34*Sterms[[2]]
    C <- 71.595
    D <- -0.077993 + 0.0028426*Sterms[[2]]
    E <- 0
    
    Ksp_calcite <- 10^(logK(A, B, C, D, E, T)) * exp(deltaPlnK(T, p, DeltaPcoeffs$Ksp_calcite))  
    attr(Ksp_calcite, "unit") <- "(mol/kg-soln)^2"     

    return(Ksp_calcite)
  }

Ksp_aragonite <- function(S, t, p=0) 
  {
    T <- T(t)
    Sterms <- Sterms(S)
    
    A <- -171.945  -  0.068393*Sterms[[2]] - 0.10018*S + 0.0059415*Sterms[[3]]
    B <- 2903.293  +    88.135*Sterms[[2]]
    C <- 71.595
    D <- -0.077993 + 0.0017276*Sterms[[2]]
    E <- 0
    
    Ksp_aragonite <- 10^(logK(A, B, C, D, E, T)) * exp(deltaPlnK(T, p, DeltaPcoeffs$Ksp_aragonite))
    attr(Ksp_aragonite, "unit") <- "(mol/kg-soln)^2"     

    return(Ksp_aragonite)
  }
