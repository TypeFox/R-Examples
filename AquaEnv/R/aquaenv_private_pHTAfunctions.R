# concentrations of species of an univalent acid
HAuni <- function(Sum,K,H)
  {
    H/(H + K)*Sum
  }

Auni <- function(Sum,K,H)
  {
    K/(H + K)*Sum
  }


# concentrations of species of a bivalent acid
H2Abi <- function(Sum, K1, K2, H)
  {
    (H^2/(H^2+K1*H+K1*K2))*Sum
  }

HAbi <- function(Sum, K1, K2, H)
  {
    (K1*H/(H^2+K1*H+K1*K2))*Sum
  }

Abi <- function(Sum, K1, K2, H)
  {
    (K1*K2/(H^2+K1*H+K1*K2))*Sum
  }


# concentrations of species of a trivalent acid
H3Atri <- function(Sum, K1, K2, K3, H)
  {
    (H^3/(H^3 + H^2*K1 + H*K1*K2 + K1*K2*K3))*Sum
  }

H2Atri <- function(Sum, K1, K2, K3, H)
  {
    (K1*H^2/(H^3 + H^2*K1 + H*K1*K2 + K1*K2*K3))*Sum
  }

HAtri <- function(Sum, K1, K2, K3, H)
  {
    (K1*K2*H/(H^3 + H^2*K1 + H*K1*K2 + K1*K2*K3))*Sum
  }

Atri <- function(Sum, K1, K2, K3, H)
  {
    (K1*K2*K3/(H^3 + H^2*K1 + H*K1*K2 + K1*K2*K3))*Sum
  }


# PRIVATE function:
# calculates [TA] from an object of class aquanenv and a given [H+]
calcTA <- function(aquaenv,                   # object of class aquaenv
                   H)                         # given [H+] in mol/kg-solution
{
  with (aquaenv,
        {
          HCO3   <- HAbi  (SumCO2,   K_CO2,   K_HCO3,           H)
          CO3    <- Abi   (SumCO2,   K_CO2,   K_HCO3,           H)

          return(HCO3 + 2*CO3 + calcTAMinor(aquaenv, H)) # the calculated [TA]
        })
}


# PRIVATE function:
# calculates minor contributions to [TA] from an object of class aquanenv and a given [H+]
calcTAMinor <- function(aquaenv,              # object of class aquaenv
                        H)                    # given [H+] in mol/kg-solution
  {
    with (aquaenv,
          {
            BOH4   <- Auni  (SumBOH3,  K_BOH3,                    H)
            OH     <- K_W  / H
            HPO4   <- HAtri (SumH3PO4, K_H3PO4, K_H2PO4,  K_HPO4, H)
            PO4    <- Atri  (SumH3PO4, K_H3PO4, K_H2PO4,  K_HPO4, H)
            SiOOH3 <- HAbi  (SumSiOH4, K_SiOH4, K_SiOOH3,         H) 
            HS     <- HAbi  (SumH2S,   K_H2S,   K_HS,             H)
            S2min  <- Abi   (SumH2S,   K_H2S,   K_HS,             H)
            NH3    <- Auni  (SumNH4,   K_NH4,                     H)
            HSO4   <- HAbi  (SumH2SO4, K_H2SO4, K_HSO4,           H)
            HF     <- HAuni (SumHF,    K_HF,                      H)
            H3PO4  <- H3Atri(SumH3PO4, K_H3PO4, K_H2PO4, K_HPO4,  H)

            return(BOH4 + OH + HPO4 + 2*PO4 + SiOOH3 + HS + 2*S2min + NH3 - H - HSO4 - HF - H3PO4)  # calculated minor contributions to [TA]
          })
  }


# PRIVATE function:
# calculates [H+]  from an object of class aquanenv and a given [TA]: first according to Follows2006, if no solution is found after Technicals$maxiter iterations, uniroot is applied
calcH_TA <- function(aquaenv,                 # object of class aquaenv     
                     TA)                      # given [TA] in mol/kg-solution
  {
    H <- c()
    aquaenv[["TA"]] <- TA
    for (z in 1:length(aquaenv))
      {
        aquaenvtemp <- as.list(as.data.frame(aquaenv)[z,])
        with (aquaenvtemp,
              {              
                Htemp <- Technicals$Hstart;  i <- 1
                while ((abs(calcTA(aquaenvtemp, Htemp) - aquaenvtemp$TA) > Technicals$Haccur) && (i <= Technicals$maxiter))
                  {
                    a <- TA - calcTAMinor(aquaenvtemp, Htemp)
                    b <- K_CO2*(a-SumCO2)
                    c <- K_CO2*K_HCO3*(a-2*SumCO2)
                   
                    Htemp <- (-b + sqrt(b^2 - (4*a*c)))/(2*a); i <- i + 1             
                    if (Htemp<0) {break}                   
                  }
                if ((Htemp < 0) || (i > Technicals$maxiter))
                  {
                    Htemp <- uniroot(function(x){calcTA(aquaenvtemp, x) - aquaenvtemp$TA},
                                     Technicals$unirootinterval, tol=Technicals$uniroottol, maxiter=Technicals$maxiter)$root
                  }
                H <<- c(H, Htemp)
              })
      }
    return(H)                                 # calculated [H+] in mol/kg-solution
  }


# PRIVATE function:
# calculates [H+]  from an object of class aquanenv and a given [CO2]: by analytically solving the resulting quadratic equation
calcH_CO2 <- function(aquaenv,                # object of class aquaenv  
                      CO2)                    # given [CO2] in mol/kg-solution
  {
    H <- c()
    aquaenv[["CO2"]] <- CO2
    for (x in 1:length(aquaenv))
      {
        aquaenvtemp <- as.list(as.data.frame(aquaenv)[x,]) 
        with (aquaenvtemp,
                {
                  a <- CO2 - SumCO2
                  b <- K_CO2*CO2
                  c <- K_CO2*K_HCO3*CO2
                  
                  H <<- c(H, ((-b-sqrt(b^2 - 4*a*c))/(2*a)))
                })
      }
    return(H)                                 # calculated [H+] in mol/kg-solution
  }
    

# PRIVATE function:
# calculates [SumCO2]  from an object of class aquanenv, a given pH, and a given [CO2]: by analytically solving the resulting equation
calcSumCO2_pH_CO2 <- function(aquaenv,        # object of class aquaenv     
                              pH,             # given pH on the free proton scale
                              CO2)            # given [CO2] in mol/kg-solution
  {
    SumCO2 <- c()
    aquaenv[["pH"]]  <- pH
    aquaenv[["CO2"]] <- CO2
    for (x in 1:length(aquaenv))
      {
        aquaenvtemp <- as.list(as.data.frame(aquaenv)[x,]) 
        with (aquaenvtemp,
                {
                  H     <- 10^{-pH}
                  denom <- H^2/(H^2 + H*K_CO2 + K_CO2*K_HCO3)
                  
                  SumCO2 <<- c(SumCO2, (CO2/denom))
                })
      }
    return(SumCO2)                            # calculated [SumCO2] in mol/kg-solution
  }


# PRIVATE function:
# calculates [SumCO2]  from an object of class aquanenv, a given pH, and a given [TA]: by analytically solving the resulting quadratic equation
calcSumCO2_pH_TA <- function(aquaenv,         # object of class aquaenv   
                             pH,              # given pH on the free proton scale
                             TA)              # given [TA] in mol/kg-solution
  {
   SumCO2 <- c()
   aquaenv[["pH"]] <- pH
   aquaenv[["TA"]] <- TA
   for (x in 1:length(aquaenv))
     {
       aquaenvtemp <- as.list(as.data.frame(aquaenv)[x,]) 
       with (aquaenvtemp,
             {
               H     <- 10^{-pH}
               c2    <- (H*K_CO2)     /(H^2 + H*K_CO2 + K_CO2*K_HCO3)
               c3    <- (K_CO2*K_HCO3)/(H^2 + H*K_CO2 + K_CO2*K_HCO3)
               numer <- TA - calcTAMinor(aquaenvtemp, H)
               denom <- c2 + 2*c3
               
               SumCO2 <<- c(SumCO2, (numer/denom))
             })
     }
   return(SumCO2)                             # calculated [SumCO2] in mol/kg-solution
 }


# PRIVATE function:
# calculates [SumCO2] from an object of class aquanenv, a given [TA], and a given [CO2]: by analytically solving the resulting quadratic equation
calcSumCO2_TA_CO2 <- function(aquaenv,        # object of class aquaenv
                              TA,             # given [TA] in mol/kg-solution
                              CO2)            # given [CO2] in mol/kg-solution
  {
    SumCO2 <- c()
    aquaenv[["TA"]] <- TA
    aquaenv[["CO2"]] <- CO2
    for (i in 1:length(aquaenv))
      {
        aquaenvtemp <- as.list(as.data.frame(aquaenv)[i,]) 
        with (aquaenvtemp,
              {
                f <- function(x)
                  {
                    c1    <- (x^2)         /(x^2 + x*K_CO2 + K_CO2*K_HCO3)
                    c2    <- (x*K_CO2)     /(x^2 + x*K_CO2 + K_CO2*K_HCO3)
                    c3    <- (K_CO2*K_HCO3)/(x^2 + x*K_CO2 + K_CO2*K_HCO3)
                             
                    return(TA - ((c2+2*c3)*(CO2/c1) + calcTAMinor(aquaenvtemp, x)))
                  }
                H     <- uniroot(f, Technicals$unirootinterval, tol=Technicals$uniroottol, maxiter=Technicals$maxiter)$root

                attr(aquaenvtemp, "class") <- "aquaenv"
                SumCO2 <<- c(SumCO2, calcSumCO2_pH_TA(aquaenvtemp, -log10(H), TA))
              })
      }
    return(SumCO2)                            # calculated [SumCO2] in mol/kg-solution
  }
