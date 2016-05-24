###########################################################################
# Preparing variables for the constructor of the class aquaenv: conversions
###########################################################################

convert <- function(x, ...)
  {
    if ((!is.null(attr(x, "class"))) && (attr(x, "class") == "aquaenv"))
      {
        return(convert.aquaenv(x, ...))
      }
    else
      {
        return(convert.standard(x, ...))
      }
  }


#########################################################################
# The class aquaenv: constructor, plotting and conversion to a data frame
#########################################################################

aquaenv <- function(S, t, p=pmax((P-Pa), gauge_p(d, lat, Pa)), P=Pa, Pa=1.01325, d=0, lat=0,
                    SumCO2=0, SumNH4=0, SumH2S=0, SumH3PO4=0, SumSiOH4=0, SumHNO3=0, SumHNO2=0, 
                    SumBOH3=NULL, SumH2SO4=NULL, SumHF=NULL,
                    TA=NULL, pH=NULL, fCO2=NULL, CO2=NULL,
                    speciation=TRUE,
                    dsa=FALSE,
                    ae=NULL,
                    from.data.frame=FALSE,
                    SumH2SO4_Koffset=0,
                    SumHF_Koffset=0,
                    revelle=FALSE,
                    skeleton=FALSE,
                    k_w=NULL,
                    k_co2=NULL,
                    k_hco3=NULL,
                    k_boh3=NULL,
                    k_hso4=NULL,
                    k_hf=NULL,
                    k1k2="roy",
                    khf="dickson",
                    khso4="dickson",
                    fCO2atm=0.000390,
                    fO2atm=0.20946)
  { 
    if (from.data.frame)
      {
        return (from.data.frame(ae))
      }
    else if (!is.null(ae))
      {
        return (cloneaquaenv(ae, TA=TA, pH=pH, k_co2=k_co2, k1k2=k1k2, khf=khf, khso4=khso4))
      }
    else
      {
        aquaenv <- list()
        attr(aquaenv, "class") <- "aquaenv"

        aquaenv[["S"]]           <- S                   ; attr(aquaenv[["S"]], "unit")           <- "p2su"        
        aquaenv[["t"]]           <- t                   ; attr(aquaenv[["t"]], "unit")           <- "deg C"
        aquaenv[["p"]]           <- p                   ; attr(aquaenv[["p"]], "unit")           <- "bar"

        aquaenv[["T"]]           <- T(t)                ; attr(aquaenv[["T"]], "unit")           <- "deg K"
                
        aquaenv[["Cl"]]          <- Cl(S)               ; attr(aquaenv[["Cl"]], "unit")          <- "permil"
        aquaenv[["I"]]           <- I(S)                ; attr(aquaenv[["I"]], "unit")           <- "mol/kg-H2O"
        
        aquaenv[["P"]]           <- p+Pa                ; attr(aquaenv[["P"]], "unit")           <- "bar"
        aquaenv[["Pa"]]          <- Pa                  ; attr(aquaenv[["Pa"]], "unit")           <- "bar"
        aquaenv[["d"]]           <- round(watdepth(P=p+Pa, lat=lat),2) ; attr(aquaenv[["d"]], "unit")      <- "m"

        aquaenv[["density"]]     <- seadensity(S=S, t=t)     ; attr(aquaenv[["density"]], "unit")      <- "kg/m3"
        
        aquaenv[["SumCO2"]]      <- SumCO2              
        aquaenv[["SumNH4"]]      <- SumNH4              ; attr(aquaenv[["SumNH4"]], "unit")       <- "mol/kg-soln"
        aquaenv[["SumH2S"]]      <- SumH2S              ; attr(aquaenv[["SumH2S"]], "unit")       <- "mol/kg-soln"
        aquaenv[["SumHNO3"]]     <- SumHNO3             ; attr(aquaenv[["SumHNO3"]], "unit")      <- "mol/kg-soln"
        aquaenv[["SumHNO2"]]     <- SumHNO2             ; attr(aquaenv[["SumHNO2"]], "unit")      <- "mol/kg-soln"
        aquaenv[["SumH3PO4"]]    <- SumH3PO4            ; attr(aquaenv[["SumH3PO4"]], "unit")     <- "mol/kg-soln"
        aquaenv[["SumSiOH4"]]    <- SumSiOH4            ; attr(aquaenv[["SumSiOH4"]], "unit")     <- "mol/kg-soln"
        if (is.null(SumBOH3))
          {
            SumBOH3 = seaconc("B", S)
          }
        else if (length(S)>1)
          {
            SumBOH3 = rep(SumBOH3, length(S))
          }
        aquaenv[["SumBOH3"]]     <- SumBOH3             ; attr(aquaenv[["SumBOH3"]], "unit")      <- "mol/kg-soln"
        if (is.null(SumH2SO4))
          {
            SumH2SO4 = seaconc("SO4", S)
          }
        else if (length(S)>1)
          {
            SumH2SO4 = rep(SumH2SO4, length(S))
          }
        aquaenv[["SumH2SO4"]]    <- SumH2SO4            ; attr(aquaenv[["SumH2SO4"]], "unit")     <- "mol/kg-soln"
        if (is.null(SumHF))
          {
            SumHF = seaconc("F", S)
          }
        else if (length(S)>1)
          {
            SumHF = rep(SumHF, length(S))
          }
        aquaenv[["SumHF"]]       <- SumHF               ;  attr(aquaenv[["SumHF"]], "unit")       <- "mol/kg-soln"

        if(!skeleton)
          {
            aquaenv[["Br"]]          <- seaconc("Br", S)    ; attr(aquaenv[["Br"]], "unit")        <- "mol/kg-soln"
            
            aquaenv[["ClConc"]]      <- seaconc("Cl", S)    ; attr(aquaenv[["ClConc"]], "unit")       <- "mol/kg-soln"
            aquaenv[["Na"]]          <- seaconc("Na", S)    ; attr(aquaenv[["Na"]], "unit")           <- "mol/kg-soln"
            aquaenv[["Mg"]]          <- seaconc("Mg", S)    ; attr(aquaenv[["Mg"]], "unit")           <- "mol/kg-soln"
            
            aquaenv[["Ca"]]          <- seaconc("Ca", S)    ; attr(aquaenv[["Ca"]], "unit")           <- "mol/kg-soln"
            
            aquaenv[["K"]]           <- seaconc("K", S)     ; attr(aquaenv[["K"]], "unit")            <- "mol/kg-soln"
            aquaenv[["Sr"]]          <- seaconc("Sr", S)    ; attr(aquaenv[["Sr"]], "unit")           <- "mol/kg-soln"
            
            
            aquaenv[["molal2molin"]] <- molal2molin(S)      ; attr(aquaenv[["molal2molin"]], "unit")  <- "(mol/kg-soln)/(mol/kg-H2O)"
            
            scaleconvs <- scaleconvert(S, t, p, SumH2SO4, SumHF, khf=khf, khso4=khso4)
            aquaenv[["free2tot"]]    <- scaleconvs$free2tot ; attr(aquaenv[["free2tot"]], "pH scale") <- "free -> tot"
            aquaenv[["free2sws"]]    <- scaleconvs$free2sws ; attr(aquaenv[["free2sws"]], "pH scale") <- "free -> sws"
            aquaenv[["tot2free"]]    <- scaleconvs$tot2free ; attr(aquaenv[["tot2free"]], "pH scale") <- "total -> free"
            aquaenv[["tot2sws"]]     <- scaleconvs$tot2sws  ; attr(aquaenv[["tot2sws"]], "pH scale")  <- "total -> sws"
            aquaenv[["sws2free"]]    <- scaleconvs$sws2free ; attr(aquaenv[["sws2free"]], "pH scale") <- "sws -> free"
            aquaenv[["sws2tot"]]     <- scaleconvs$sws2tot  ; attr(aquaenv[["sws2tot"]], "pH scale")  <- "sws -> total"
            
        
            aquaenv[["K0_CO2"]]      <- K0_CO2  (S, t)
            aquaenv[["K0_O2"]]       <- K0_O2   (S, t)

            aquaenv[["fCO2atm"]]     <- fCO2atm             ; attr(aquaenv[["fCO2atm"]], "unit")      <- "atm"
            aquaenv[["fO2atm"]]      <- fO2atm              ; attr(aquaenv[["fO2atm"]], "unit")       <- "atm"
            
            aquaenv[["CO2_sat"]]     <- aquaenv[["fCO2atm"]] * aquaenv[["K0_CO2"]] ; attr(aquaenv[["CO2_sat"]], "unit")  <- "mol/kg-soln"
            aquaenv[["O2_sat"]]      <- aquaenv[["fO2atm"]]  * aquaenv[["K0_O2"]]  ; attr(aquaenv[["O2_sat"]],  "unit")  <- "mol/kg-soln"
          }

        len <- max(length(t), length(S), length(p))
        
        if (is.null(k_w))
          {
            aquaenv[["K_W"]]         <- K_W(S, t, p, SumH2SO4 + SumH2SO4_Koffset, SumHF + SumHF_Koffset, khf=khf, khso4=khso4)
          }
        else
          {
            aquaenv[["K_W"]]         <- eval(att(rep(k_w,len)))
          }

        if (is.null(k_hso4))
          {    
            aquaenv[["K_HSO4"]]      <- K_HSO4(S, t, p, khso4)
          }
        else
          {
            aquaenv[["K_HSO4"]]      <- eval(att(rep(k_hso4,len)))
          }

        if (is.null(k_hf))
          {
            aquaenv[["K_HF"]]    <- K_HF(S, t, p, SumH2SO4=(SumH2SO4 + SumH2SO4_Koffset), SumHF=(SumHF + SumHF_Koffset), khf=khf, khso4=khso4)
          }
        else
        {
          aquaenv[["K_HF"]]          <- eval(att(rep(k_hf,len)))
        }

        if (is.null(k_co2))
          {
            aquaenv[["K_CO2"]]    <- K_CO2 (S, t, p, SumH2SO4 + SumH2SO4_Koffset, SumHF + SumHF_Koffset, k1k2=k1k2, khf=khf, khso4=khso4)
          }
        else
          {
            aquaenv[["K_CO2"]]       <- eval(att(rep(k_co2,len)))
          }

        if (is.null(k_hco3))
          {
            aquaenv[["K_HCO3"]]  <- K_HCO3(S, t, p, SumH2SO4 + SumH2SO4_Koffset, SumHF + SumHF_Koffset, k1k2=k1k2, khf=khf, khso4=khso4)
          }
        else
          {
            aquaenv[["K_HCO3"]]      <- eval(att(rep(k_hco3,len)))
          }

        if (is.null(k_boh3))
          {
            aquaenv[["K_BOH3"]]      <- K_BOH3(S, t, p, SumH2SO4 + SumH2SO4_Koffset, SumHF + SumHF_Koffset, khf=khf, khso4=khso4)
          }
        else
          {
            aquaenv[["K_BOH3"]]      <- eval(att(rep(k_boh3,len)))
          }
        
        aquaenv[["K_NH4"]]       <- K_NH4   (S, t, p, SumH2SO4 + SumH2SO4_Koffset, SumHF + SumHF_Koffset, khf=khf, khso4=khso4)
        aquaenv[["K_H2S"]]       <- K_H2S   (S, t, p, SumH2SO4 + SumH2SO4_Koffset, SumHF + SumHF_Koffset, khf=khf, khso4=khso4)  
        aquaenv[["K_H3PO4"]]     <- K_H3PO4 (S, t, p, SumH2SO4 + SumH2SO4_Koffset, SumHF + SumHF_Koffset, khf=khf, khso4=khso4)
        aquaenv[["K_H2PO4"]]     <- K_H2PO4 (S, t, p, SumH2SO4 + SumH2SO4_Koffset, SumHF + SumHF_Koffset, khf=khf, khso4=khso4)
        aquaenv[["K_HPO4"]]      <- K_HPO4  (S, t, p, SumH2SO4 + SumH2SO4_Koffset, SumHF + SumHF_Koffset, khf=khf, khso4=khso4)
        aquaenv[["K_SiOH4"]]     <- K_SiOH4 (S, t, p, SumH2SO4 + SumH2SO4_Koffset, SumHF + SumHF_Koffset, khf=khf, khso4=khso4)
        aquaenv[["K_SiOOH3"]]    <- K_SiOOH3(S, t, p, SumH2SO4 + SumH2SO4_Koffset, SumHF + SumHF_Koffset, khf=khf, khso4=khso4)  
        
        aquaenv[["K_HNO2"]]      <- PhysChemConst$K_HNO2    ; attr(aquaenv[["K_HNO2"]], "unit")       <- "mol/kg-soln; mol/kg-H2O; mol/l"
        aquaenv[["K_HNO3"]]      <- PhysChemConst$K_HNO3    ; attr(aquaenv[["K_HNO3"]], "unit")       <- "mol/kg-soln; mol/kg-H2O; mol/l"
        aquaenv[["K_H2SO4"]]     <- PhysChemConst$K_H2SO4   ; attr(aquaenv[["K_H2SO4"]], "unit")      <- "mol/kg-soln; mol/kg-H2O; mol/l"
        aquaenv[["K_HS"]]        <- PhysChemConst$K_HS      ; attr(aquaenv[["K_HS"]], "unit")         <- "mol/kg-soln; mol/kg-H2O; mol/l"

        if (!skeleton)
          {
            aquaenv[["Ksp_calcite"]]    <- Ksp_calcite(S, t, p)
            aquaenv[["Ksp_aragonite"]]  <- Ksp_aragonite(S, t, p)
          }
               
        if (!(is.null(TA) && is.null(pH) && is.null(fCO2) && is.null(CO2)))
          {
            if (is.null(SumCO2))
              {
                if (!is.null(pH))
                  {
                    if (!is.null(fCO2))
                      {
                        CO2    <- aquaenv[["K0_CO2"]] * fCO2
                        SumCO2 <- calcSumCO2_pH_CO2(aquaenv, pH, CO2)
                        TA     <- calcTA(c(aquaenv, list(SumCO2, "SumCO2")), 10^(-pH))
                      }
                    else if (!is.null(CO2))
                      {
                        SumCO2 <- calcSumCO2_pH_CO2(aquaenv, pH, CO2)
                        fCO2   <- CO2 / aquaenv[["K0_CO2"]]
                        TA     <- calcTA(c(aquaenv, list(SumCO2, "SumCO2")), 10^(-pH))
                      }
                    else if (!is.null(TA))
                      {
                        SumCO2 <- calcSumCO2_pH_TA(aquaenv, pH, TA)
                        CO2    <- H2Abi(SumCO2, aquaenv[["K_CO2"]], aquaenv[["K_HCO3"]], 10^(-pH))
                        fCO2   <- CO2 / aquaenv[["K0_CO2"]]
                      }
                    else
                      {
                        print(paste("Error! Underdetermined system!"))
                        print("To calculate [SumCO2], please enter one pair of: (pH/CO2), (pH,fCO2), (TA,pH), (TA/CO2), (TA/fCO2)")
                        return(NULL) 
                      }
                  }
                else if (!is.null(TA))
                  {
                    if (!is.null(fCO2))
                      {
                        CO2    <- aquaenv[["K0_CO2"]] * fCO2
                        SumCO2 <- calcSumCO2_TA_CO2(aquaenv, TA, CO2) 
                        pH     <- -log10(calcH_TA(c(aquaenv, list(SumCO2, "SumCO2")), TA))
                      }
                    else if (!is.null(CO2))
                      {
                        SumCO2 <- calcSumCO2_TA_CO2(aquaenv, TA, CO2)
                        fCO2   <- CO2 / aquaenv[["K0_CO2"]]
                        pH     <- -log10(calcH_TA(c(aquaenv, list(SumCO2, "SumCO2")), TA))
                      }
                    else
                      {
                        print(paste("Error! Underdetermined system!"))
                        print("To calculate [SumCO2], please enter one pair of: (pH/CO2), (pH,fCO2), (TA,pH), (TA/CO2), (TA/fCO2)")
                        return(NULL) 
                      }
                  }
                else
                  {
                    print(paste("Error! Underdetermined system!"))
                    print("To calculate [SumCO2], please enter one pair of: (pH/CO2), (pH,fCO2), (TA,pH), (TA/CO2), (TA/fCO2)")
                    return(NULL) 
                  }
              }
            else
              {        
                if (is.null(fCO2) && is.null(CO2) && !is.null(pH) && !is.null(TA))
                  {
                    TAcalc <- calcTA(aquaenv, 10^(-pH))
                    print(paste("Error! Overdetermined system: entered TA:", TA, ", calculated TA:", TAcalc))
                    print("Please enter only one of: pH, TA, CO2, or fCO2.")
                    return(NULL)
                  }
                if ((!is.null(pH) || !is.null(TA)) && (is.null(fCO2) && is.null(CO2)))
                  {
                    if (!is.null(pH))
                      {
                        H              <- 10^(-pH)
                        TA             <- calcTA(aquaenv, H)
                        attributes(TA) <- NULL
                      }
                    else 
                      {
                        pH             <- -log10(calcH_TA(aquaenv, TA))
                      }
                    CO2  <- H2Abi(SumCO2, aquaenv[["K_CO2"]], aquaenv[["K_HCO3"]], 10^(-pH))
                    attributes(CO2) <- NULL
                    fCO2 <- CO2 / aquaenv[["K0_CO2"]]
                    attributes(fCO2) <- NULL
                  }
                else if (!is.null(fCO2) && !is.null(CO2) && is.null(pH) && is.null(TA))
                  {
                    fCO2calc <- CO2 / aquaenv[["K0_CO2"]]
                    print(paste("Error! Overdetermined system: entered fCO2:", fCO2, ", calculated fCO2:", fCO2calc))
                    print("Please enter only one of: pH, TA, CO2, or fCO2.")
                    return(NULL)
                  }
                else if ((!is.null(fCO2) || !is.null(CO2)) && (is.null(pH) && is.null(TA)))
                  {
                    if(!is.null(fCO2))
                      {
                        CO2 <- aquaenv[["K0_CO2"]] * fCO2
                        attributes(CO2) <- NULL
                      }
                    else
                      {
                        fCO2 <- CO2 / aquaenv[["K0_CO2"]]
                        attributes(fCO2) <- NULL
                      }
                    H   <- calcH_CO2(aquaenv, CO2)
                    TA  <- calcTA(aquaenv, H)
                    pH <- -log10(H)
                  }
                else if ((!is.null(fCO2) || !is.null(CO2)) && !is.null(pH) && is.null(TA))
                  {
                    if (!is.null(fCO2))
                      {
                        CO2 <- aquaenv[["K0_CO2"]] * fCO2
                      }
                    pHcalc <- -log10(calcH_CO2(aquaenv, CO2))
                    print(paste("Error! Overdetermined system: entered pH:", pH, ", calculated pH:", pHcalc))
                    print("Please enter only one of: pH, TA, CO2, or fCO2.")
                    return(NULL)
                  }
                else if ((!is.null(fCO2) || !is.null(CO2)) && is.null(pH) && !is.null(TA))
                  {
                    if (!is.null(fCO2))
                      {
                        CO2 <- aquaenv[["K0_CO2"]] * fCO2
                      }
                    H   <- calcH_CO2(aquaenv, CO2)
                    TAcalc <- calcTA(aquaenv, H)
                    print(paste("Error! Overdetermined system: entered TA:", TA, ", calculated TA:", TAcalc))
                    print("Please enter only one of: pH, TA, CO2, or fCO2.")
                    return(NULL)
                  }
              }
            
            aquaenv[["TA"]]     <- TA     ; attr(aquaenv[["TA"]], "unit")      <- "mol/kg-soln"
            aquaenv[["pH"]]     <- pH     ; attr(aquaenv[["pH"]], "pH scale")  <- "free"
            aquaenv[["SumCO2"]] <- SumCO2 ; attr(aquaenv[["SumCO2"]], "unit")  <- "mol/kg-soln"
            aquaenv[["fCO2"]]   <- fCO2   ; attr(aquaenv[["fCO2"]], "unit")    <- "atm"
            aquaenv[["CO2"]]    <- CO2    ; attr(aquaenv[["CO2"]], "unit")     <- "mol/kg-soln"
            
            
            if (speciation)
              {
                H <- 10^(-pH)
                
                with (aquaenv,
                      {
                        aquaenv[["HCO3"]]   <<- HAbi  (SumCO2,   K_CO2,   K_HCO3,          H) ; attr(aquaenv[["HCO3"]], "unit")    <<- "mol/kg-soln"
                        aquaenv[["CO3"]]    <<- Abi   (SumCO2,   K_CO2,   K_HCO3,          H) ; attr(aquaenv[["CO3"]], "unit")     <<- "mol/kg-soln"
                        
                        aquaenv[["BOH3"]]   <<- HAuni (SumBOH3,  K_BOH3,                   H) ; attr(aquaenv[["BOH3"]], "unit")    <<- "mol/kg-soln"
                        aquaenv[["BOH4"]]   <<- Auni  (SumBOH3,  K_BOH3,                   H) ; attr(aquaenv[["BOH4"]], "unit")    <<- "mol/kg-soln"
                        
                        aquaenv[["OH"]]     <<- K_W  / H                                      ; attr(aquaenv[["OH"]], "unit")      <<- "mol/kg-soln"
                        
                        aquaenv[["H3PO4"]]  <<- H3Atri(SumH3PO4, K_H3PO4, K_H2PO4, K_HPO4, H) ; attr(aquaenv[["H3PO4"]], "unit")   <<- "mol/kg-soln"
                        aquaenv[["H2PO4"]]  <<- H2Atri(SumH3PO4, K_H3PO4, K_H2PO4, K_HPO4, H) ; attr(aquaenv[["H2PO4"]], "unit")   <<- "mol/kg-soln"
                        aquaenv[["HPO4"]]   <<- HAtri (SumH3PO4, K_H3PO4, K_H2PO4, K_HPO4, H) ; attr(aquaenv[["HPO4"]], "unit")    <<- "mol/kg-soln"
                        aquaenv[["PO4"]]    <<- Atri  (SumH3PO4, K_H3PO4, K_H2PO4, K_HPO4, H) ; attr(aquaenv[["PO4"]], "unit")     <<- "mol/kg-soln"
                        
                        aquaenv[["SiOH4"]]  <<- H2Abi (SumSiOH4, K_SiOH4, K_SiOOH3,        H) ; attr(aquaenv[["SiOH4"]], "unit")   <<- "mol/kg-soln"
                        aquaenv[["SiOOH3"]] <<- HAbi  (SumSiOH4, K_SiOH4, K_SiOOH3,        H) ; attr(aquaenv[["SiOOH3"]], "unit")  <<- "mol/kg-soln"
                        aquaenv[["SiO2OH2"]]<<-  Abi  (SumSiOH4, K_SiOH4, K_SiOOH3,        H) ; attr(aquaenv[["SiO2OH2"]], "unit") <<- "mol/kg-soln"
                        
                        aquaenv[["H2S"]]    <<- H2Abi (SumH2S,   K_H2S,   K_HS,            H) ; attr(aquaenv[["H2S"]], "unit")     <<- "mol/kg-soln"
                        aquaenv[["HS"]]     <<- HAbi  (SumH2S,   K_H2S,   K_HS,            H) ; attr(aquaenv[["HS"]], "unit")      <<- "mol/kg-soln"
                        aquaenv[["S2min"]]  <<- Abi   (SumH2S,   K_H2S,   K_HS,            H) ; attr(aquaenv[["S2min"]], "unit")   <<- "mol/kg-soln"
                        
                        aquaenv[["NH4"]]    <<- HAuni (SumNH4,   K_NH4,                    H) ; attr(aquaenv[["NH4"]], "unit")     <<- "mol/kg-soln"
                        aquaenv[["NH3"]]    <<- Auni  (SumNH4,   K_NH4,                    H) ; attr(aquaenv[["NH3"]], "unit")     <<- "mol/kg-soln"
                        
                        aquaenv[["H2SO4"]]  <<- H2Abi (SumH2SO4, K_H2SO4, K_HSO4,          H) ; attr(aquaenv[["H2SO4"]], "unit")   <<- "mol/kg-soln"
                        aquaenv[["HSO4"]]   <<- HAbi  (SumH2SO4, K_H2SO4, K_HSO4,          H) ; attr(aquaenv[["HSO4"]], "unit")    <<- "mol/kg-soln"
                        aquaenv[["SO4"]]    <<- Abi   (SumH2SO4, K_H2SO4, K_HSO4,          H) ; attr(aquaenv[["SO4"]], "unit")     <<- "mol/kg-soln"
                        
                        aquaenv[["HF"]]     <<- HAuni (SumHF,    K_HF,                     H) ; attr(aquaenv[["HF"]], "unit")      <<- "mol/kg-soln"
                        aquaenv[["F"]]      <<- Auni  (SumHF,    K_HF,                     H) ; attr(aquaenv[["F"]], "unit")       <<- "mol/kg-soln"

                        aquaenv[["HNO3"]]   <<- HAuni (SumHNO3,  K_HNO3,                   H) ; attr(aquaenv[["HNO3"]], "unit")    <<- "mol/kg-soln"
                        aquaenv[["NO3"]]    <<- Auni  (SumHNO3,  K_HNO3,                   H) ; attr(aquaenv[["NO3"]], "unit")     <<- "mol/kg-soln"
                        
                        aquaenv[["HNO2"]]   <<- HAuni (SumHNO2,  K_HNO2,                   H) ; attr(aquaenv[["HNO2"]], "unit")    <<- "mol/kg-soln"
                        aquaenv[["NO2"]]    <<- Auni  (SumHNO2,  K_HNO2,                   H) ; attr(aquaenv[["NO2"]], "unit")     <<- "mol/kg-soln"


                        aquaenv[["omega_calcite"]]               <<- aquaenv[["Ca"]]*aquaenv[["CO3"]]/aquaenv[["Ksp_calcite"]]
                        attributes(aquaenv[["omega_calcite"]])   <<- NULL
                        aquaenv[["omega_aragonite"]]             <<- aquaenv[["Ca"]]*aquaenv[["CO3"]]/aquaenv[["Ksp_aragonite"]]
                        attributes(aquaenv[["omega_aragonite"]]) <<- NULL
                      })

                if (revelle)
                  {
                    aquaenv[["revelle"]]             <- revelle(aquaenv)
                    attributes(aquaenv[["revelle"]]) <- NULL
                  }
              }
            if (dsa)
              {
                with (aquaenv,
                      {
                        if (!((length(SumCO2)==1) && (SumCO2==0)))
                          {
                            aquaenv[["c1"]]  <<- CO2 /SumCO2      ; attributes(aquaenv[["c1"]])  <<- NULL
                            aquaenv[["c2"]]  <<- HCO3/SumCO2      ; attributes(aquaenv[["c2"]])  <<- NULL
                            aquaenv[["c3"]]  <<- CO3 /SumCO2      ; attributes(aquaenv[["c3"]])  <<- NULL
                            aquaenv[["dTAdSumCO2"]] <<- aquaenv[["c2"]] + 2*aquaenv[["c3"]]
                            attr(aquaenv[["dTAdSumCO2"]], "unit") <<- "(mol-TA/kg-soln)/(mol-SumCO2/kg-soln)"  
                          }
                        if(!((length(SumBOH3)==1) && (SumBOH3==0)))
                          {
                            aquaenv[["b1"]]  <<- BOH3/SumBOH3     ; attributes(aquaenv[["b1"]])  <<- NULL
                            aquaenv[["b2"]]  <<- BOH4/SumBOH3     ; attributes(aquaenv[["b2"]])  <<- NULL
                            aquaenv[["dTAdSumBOH3"]] <<- aquaenv[["b2"]]
                            attr(aquaenv[["dTAdSumBOH3"]], "unit") <<- "(mol-TA/kg-soln)/(mol-SumBOH3/kg-soln)"
                          }
                        if(!((length(SumH3PO4)==1) && (SumH3PO4==0)))
                          {
                            aquaenv[["p1"]]  <<- H3PO4/SumH3PO4   ; attributes(aquaenv[["p1"]])  <<- NULL
                            aquaenv[["p2"]]  <<- H2PO4/SumH3PO4   ; attributes(aquaenv[["p2"]])  <<- NULL
                            aquaenv[["p3"]]  <<- HPO4 /SumH3PO4   ; attributes(aquaenv[["p3"]])  <<- NULL
                            aquaenv[["p4"]]  <<- PO4  /SumH3PO4   ; attributes(aquaenv[["p4"]])  <<- NULL
                            aquaenv[["dTAdSumH3PO4"]] <<- aquaenv[["p3"]] + 2*aquaenv[["p4"]] - aquaenv[["p1"]]
                            attr(aquaenv[["dTAdSumH3PO4"]], "unit") <<- "(mol-TA/kg-soln)/(mol-SumH3PO4/kg-soln)"
                          }
                        if(!((length(SumSiOH4)==1) && (SumSiOH4==0)))
                          {
                            aquaenv[["si1"]] <<- SiOH4  /SumSiOH4 ; attributes(aquaenv[["si1"]]) <<- NULL
                            aquaenv[["si2"]] <<- SiOOH3 /SumSiOH4 ; attributes(aquaenv[["si2"]]) <<- NULL
                            aquaenv[["si3"]] <<- SiO2OH2/SumSiOH4 ; attributes(aquaenv[["si3"]]) <<- NULL
                            aquaenv[["dTAdSumSumSiOH4"]] <<- aquaenv[["si2"]]
                            attr(aquaenv[["dTAdSumSumSiOH4"]], "unit") <<- "(mol-TA/kg-soln)/(mol-SumSiOH4/kg-soln)"
                          }
                        if(!((length(SumH2S)==1) && (SumH2S==0)))
                          {
                            aquaenv[["s1"]]  <<- H2S  /SumH2S     ; attributes(aquaenv[["s1"]])  <<- NULL
                            aquaenv[["s2"]]  <<- HS   /SumH2S     ; attributes(aquaenv[["s2"]])  <<- NULL
                            aquaenv[["s3"]]  <<- S2min/SumH2S     ; attributes(aquaenv[["s3"]])  <<- NULL
                            aquaenv[["dTAdSumH2S"]] <<- aquaenv[["s2"]] + 2*aquaenv[["s3"]]
                            attr(aquaenv[["dTAdSumH2S"]], "unit") <<- "(mol-TA/kg-soln)/(mol-SumH2S/kg-soln)"
                          }
                        if(!((length(SumNH4)==1) && (SumNH4==0)))
                          {
                            aquaenv[["n1"]]  <<- NH4/SumNH4       ; attributes(aquaenv[["n1"]])  <<- NULL
                            aquaenv[["n2"]]  <<- NH3/SumNH4       ; attributes(aquaenv[["n2"]])  <<- NULL
                            aquaenv[["dTAdSumNH4"]] <<- aquaenv[["n2"]]
                            attr(aquaenv[["dTAdSumNH4"]], "unit") <<- "(mol-TA/kg-soln)/(mol-SumNH4/kg-soln)"
                          }
                        if(!((length(SumH2SO4)==1) && (SumH2SO4==0)))
                          {
                            aquaenv[["so1"]] <<- H2SO4/SumH2SO4   ; attributes(aquaenv[["so1"]])  <<- NULL
                            aquaenv[["so2"]] <<- HSO4 /SumH2SO4   ; attributes(aquaenv[["so2"]])  <<- NULL
                            aquaenv[["so3"]] <<- SO4  /SumH2SO4   ; attributes(aquaenv[["so3"]])  <<- NULL
                            aquaenv[["dTAdSumH2SO4"]] <<- -aquaenv[["so2"]]
                            attr(aquaenv[["dTAdSumH2SO4"]], "unit") <<- "(mol-TA/kg-soln)/(mol-SumH2SO4/kg-soln)"
                          }
                        if(!((length(SumHF)==1) && (SumHF==0)))
                          {
                            aquaenv[["f1"]]  <<- HF/SumHF         ; attributes(aquaenv[["f1"]])  <<- NULL
                            aquaenv[["f2"]]  <<- F /SumHF         ; attributes(aquaenv[["f2"]])  <<- NULL
                            aquaenv[["dTAdSumHF"]] <<- -aquaenv[["f1"]]
                            attr(aquaenv[["dTAdSumHF"]], "unit") <<- "(mol-TA/kg-soln)/(mol-SumHF/kg-soln)"
                          }
                        if(!((length(SumHNO3)==1) && (SumHNO3==0)))
                          {
                            aquaenv[["na1"]] <<- HNO3/SumHNO3     ; attributes(aquaenv[["na1"]])  <<- NULL
                            aquaenv[["na2"]] <<- NO3 /SumHNO3     ; attributes(aquaenv[["na2"]])  <<- NULL
                            aquaenv[["dTAdSumHNO3"]] <<- -aquaenv[["na1"]]
                            attr(aquaenv[["dTAdSumHNO3"]], "unit") <<- "(mol-TA/kg-soln)/(mol-SumHNO3/kg-soln)"
                          }
                        if(!((length(SumHNO2)==1) && (SumHNO2==0)))
                          {
                            aquaenv[["ni1"]] <<- HNO2/SumHNO2     ; attributes(aquaenv[["ni1"]])  <<- NULL
                            aquaenv[["ni2"]] <<- NO2 /SumHNO2     ; attributes(aquaenv[["ni2"]])  <<- NULL
                            aquaenv[["dTAdSumHNO2"]] <<- -aquaenv[["ni1"]]
                            attr(aquaenv[["dTAdSumHNO2"]], "unit") <<- "(mol-TA/kg-soln)/(mol-SumHNO2/kg-soln)"
                          }
                      })

                
                aquaenv[["dTAdH"]]            <- dTAdH(aquaenv)             ; attr(aquaenv[["dTAdH"]],            "unit") <- "(mol-TA/kg-soln)/(mol-H/kg-soln)"             
                aquaenv[["dTAdKdKdS"]]        <- dTAdKdKdS(aquaenv)         ; attr(aquaenv[["dTAdKdKdS"]],        "unit") <- "(mol-TA/kg-soln)/\"psu\""
                aquaenv[["dTAdKdKdT"]]        <- dTAdKdKdT(aquaenv)         ; attr(aquaenv[["dTAdKdKdT"]],        "unit") <- "(mol-TA/kg-soln)/K"
                aquaenv[["dTAdKdKdp"]]        <- dTAdKdKdp(aquaenv)         ; attr(aquaenv[["dTAdKdKdp"]],        "unit") <- "(mol-TA/kg-soln)/m"
                aquaenv[["dTAdKdKdSumH2SO4"]] <- dTAdKdKdSumH2SO4(aquaenv)  ; attr(aquaenv[["dTAdKdKdSumH2SO4"]], "unit") <- "(mol-TA/kg-soln)/(mol-SumH2SO4/kg-soln)"
                aquaenv[["dTAdKdKdSumHF"]]    <- dTAdKdKdSumHF(aquaenv)     ; attr(aquaenv[["dTAdKdKdSumHF"]],    "unit") <- "(mol-TA/kg-soln)/(mol-SumHF/kg-soln)"
              }
          }
        return(aquaenv)
      }
  }


plot.aquaenv <- function(x, xval, what=NULL, bjerrum=FALSE, cumulative=FALSE, newdevice=TRUE, setpar=TRUE, device="x11", ...)
  {
    if ((!bjerrum) && (!cumulative))
      {
        if (is.null(what))
          {
            plotall(aquaenv=x, xval=xval, newdevice=newdevice, setpar=setpar, device=device, ...)
          }
        else 
          {
            selectplot(aquaenv=x, xval=xval, what=what, newdevice=newdevice, setpar=setpar, device=device, ...)
          }
      }
    else if (bjerrum)
      {
        bjerrumplot(aquaenv=x, what=what, newdevice=newdevice, setpar=setpar,  device=device, ...)
      }
    else
      {
        cumulativeplot(aquaenv=x, xval=xval, what=what, newdevice=newdevice, setpar=setpar, device=device, ...)
      }

    if (!(device=="x11"))
      {
        dev.off()
      }
  } 


    
as.data.frame.aquaenv <- function(x, ...)
  {
    len <- length(x)
    temp  <- list()
    for (e in x)
      {
        if (length(e) > 1)
          {
            temp <- c(temp, list(e))
          }
        else
          {
            temp <- c(temp, list(rep(e,len)))
          }
      }
    names(temp) <- names(x)

    return(as.data.frame(temp))
  }



