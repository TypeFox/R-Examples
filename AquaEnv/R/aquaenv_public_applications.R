# function:
# creates an object of class aquaenv which contains a titration simulation
titration <- function(aquaenv,                # an object of type aquaenv: minimal definition, contains all information about the system: T, S, d, total concentrations of nutrients etc (Note that it is possible to give values for SumBOH4, SumHSO4, and SumHF in the sample other than the ones calculated from salinity)
                      mass_sample,            # the mass of the sample solution in kg
                      mass_titrant,           # the total mass of the added titrant solution in kg
                      conc_titrant,           # the concentration of the titrant solution in mol/kg-soln
                      S_titrant=NULL,         # the salinity of the titrant solution, if not supplied it is assumed that the titrant solution has the same salinity as the sample solution
                      steps,                  # the amount of steps the mass of titrant is added in 
                      type="HCl",             # the type of titrant: either "HCl" or "NaOH"
                      seawater_titrant=FALSE, # is the titrant based on natural seawater? (does it contain SumBOH4, SumHSO4, and SumHF in the same proportions as seawater, i.e., correlated to S?); Note that you can only assume a seawater based titrant (i.e. SumBOH4, SumHSO4, and SumHF ~ S) or a water based titrant (i.e. SumBOH4, SumHSO4, and SumHF = 0). It is not possible to give values for SumBOH4, SumHSO4, and SumHF of the titrant.
                      k_w=NULL,               # a fixed K_W can be specified
                      k_co2=NULL,             # a fixed K_CO2 can be specified: used for TA fitting: give a K_CO2 and NOT calculate it from T and S: i.e. K_CO2 can be fitted in the routine as well
                      k_hco3=NULL,            # a fixed K_HCO3 can be specified
                      k_boh3=NULL,            # a fixed K_BOH3 can be specified
                      k_hso4=NULL,            # a fixed K_HSO4 can be specified
                      k_hf=NULL,              # a fixed K_HF can be specified
                      k1k2="roy",             # either "roy" (default, Roy1993a) or "lueker" (Lueker2000, calculated with seacarb) for K\_CO2 and K\_HCO3
                      khf="dickson")          # either "dickson" (default, Dickson1979a) or "perez" (Perez1987a, calculated with seacarb) for K\_HF}
  {
    concstep <- function(conc, oldmass, newmass)
      {
        return((conc*oldmass)/(newmass))
      }

    if (is.null(S_titrant))
      {
        S_titrant <- aquaenv$S
      }
    
    directionTAchange   <- switch(type, HCl  = -1, NaOH = +1)
    TAmolchangeperstep  <- (conc_titrant*mass_titrant)/steps
    masschangeperstep   <- mass_titrant/steps

    aquaenv[["mass"]]   <- mass_sample
    
    for (i in 1:steps)
      {
        with(aquaenv,
             {
               i           <- length(mass)
               
               newmass     <- mass[[i]] + masschangeperstep
               
               SumCO2      <- concstep(SumCO2[[i]],   mass[[i]], newmass)
               SumNH4      <- concstep(SumNH4[[i]],   mass[[i]], newmass)
               SumH2S      <- concstep(SumH2S[[i]],   mass[[i]], newmass)
               SumH3PO4    <- concstep(SumH3PO4[[i]], mass[[i]], newmass)
               SumSiOH4    <- concstep(SumSiOH4[[i]], mass[[i]], newmass)
               SumHNO3     <- concstep(SumHNO3[[i]],  mass[[i]], newmass)
               SumHNO2     <- concstep(SumHNO2[[i]],  mass[[i]], newmass)
               if (seawater_titrant)
                 {
                   SumBOH3  <- (SumBOH3[[i]]*mass[[i]]  + seaconc("B", S_titrant)*masschangeperstep)/newmass
                   SumH2SO4 <- (SumH2SO4[[i]]*mass[[i]] + seaconc("SO4", S_titrant)*masschangeperstep)/newmass
                   SumHF    <- (SumHF[[i]]*mass[[i]]    + seaconc("F", S_titrant)*masschangeperstep)/newmass
                 }
               else
                 {
                   SumBOH3  <- concstep(SumBOH3[[i]],  mass[[i]], newmass)
                   SumH2SO4 <- concstep(SumH2SO4[[i]], mass[[i]], newmass)
                   SumHF    <- concstep(SumHF[[i]],    mass[[i]], newmass)
                 }
               
               S           <- (S[[i]]*mass[[i]] + S_titrant*masschangeperstep)/newmass
               
               TA          <- (TA[[i]]*mass[[i]] + (directionTAchange * TAmolchangeperstep))/newmass

               aquaenvtemp <- aquaenv(S=S, t=t[[i]], p=p[[i]], SumCO2=SumCO2, SumNH4=SumNH4, SumH2S=SumH2S, SumH3PO4=SumH3PO4, SumSiOH4=SumSiOH4, SumHNO3=SumHNO3, SumHNO2=SumHNO2,
                                      SumBOH3=SumBOH3, SumH2SO4=SumH2SO4, SumHF=SumHF, TA=TA,
                                      speciation=(!is.null(aquaenv$HCO3)), skeleton=(is.null(aquaenv$Na)), revelle=(!is.null(aquaenv$revelle)), dsa=(!is.null(aquaenv$dTAdH)),
                                      k_w=k_w, k_co2=k_co2, k_hco3=k_hco3, k_boh3=k_boh3, k_hso4=k_hso4, k_hf=k_hf, k1k2=k1k2, khf=khf)

               aquaenvtemp[["mass"]] <- newmass
               aquaenv               <<- merge(aquaenv, aquaenvtemp)              
             })
        }

    aquaenv[["delta_conc_titrant"]]   <- ((0:steps)*TAmolchangeperstep)/aquaenv[["mass"]] ; attr(aquaenv[["delta_conc_titrant"]], "unit")  <- "mol/kg-soln"
    aquaenv[["delta_mass_titrant"]]   <- (0:steps)*masschangeperstep                      ; attr(aquaenv[["delta_mass_titrant"]], "unit")  <- "kg"
    aquaenv[["delta_moles_titrant"]]  <- (0:steps)*TAmolchangeperstep                     ; attr(aquaenv[["delta_moles_titrant"]], "unit") <- "mol"

    return(aquaenv)                           # object of class aquaenv which contains a titration simulation
  }



# PUBLIC function:
# calculates [TA] and [SumCO2] from a titration curve using an optimization procedure (nls.lm from R package minpack.lm)
TAfit <- function(ae,                         # an object of type aquaenv: minimal definition, contains all information about the system: T, S, d, total concentrations of nutrients etc 
                  titcurve,                   # a table containing the titration curve: basically a series of tuples of added titrant solution mass and pH values (pH on free proton scale) or E values in V
                  conc_titrant,               # concentration of the titrant solution in mol/kg-soln
                  mass_sample,                # the mass of the sample solution in kg
                  S_titrant=NULL,             # the salinity of the titrant solution, if not supplied it is assumed that the titrant solution has the same salinity as the sample solution
                  TASumCO2guess=0.0025,       # a first guess for [TA] and [SumCO2] to be used as initial values for the optimization procedure
                  E0guess=0.4,                # first guess for E0 in V
                  type="HCl",                 # the type of titrant: either "HCl" or "NaOH"
                  Evals=FALSE,                # are the supplied datapoints pH or E (V) values?
                  electrode_polarity="pos",   # either "pos" or "neg": how is the polarity of the Electrode: E = E0 -(RT/F)ln(H^+) ("pos") or -E = E0 -(RT/F)ln(H^+) ("neg")
                  K_CO2fit=FALSE,             # should K_CO2 be fitted as well?
                  equalspaced=TRUE,           # are the mass values of titcurve equally spaced?
                  seawater_titrant=FALSE,     # is the titrant based on natural seawater? (does it contain SumBOH4, SumHSO4, and SumHF in the same proportions as seawater, i.e., correlated to S?); Note that you can only assume a seawater based titrant (i.e. SumBOH4, SumHSO4, and SumHF ~ S) or a water based titrant (i.e. SumBOH4, SumHSO4, and SumHF = 0). It is not possible to give values for SumBOH4, SumHSO4, and SumHF of the titrant.
                  pHscale="free",             # either "free", "total", "sws" or "nbs": if the titration curve contains pH data: on which scale is it measured?
                  debug=FALSE,                # debug mode: the last simulated titration tit, the converted pH profile calc, and the nls.lm output out are made global variables for investigation and plotting
                  k_w=NULL,                   # a fixed K_W can be specified
                  k_co2=NULL,                 # a fixed K_CO2 can be specified
                  k_hco3=NULL,                # a fixed K_HCO3 can be specified
                  k_boh3=NULL,                # a fixed K_BOH3 can be specified
                  k_hso4=NULL,                # a fixed K_HSO4 can be specified
                  k_hf=NULL,                  # a fixed K_HF can be specified
                  nlscontrol=nls.lm.control(),# nls.lm.control() can be specified
                  verbose=FALSE,              # verbose mode: show the traject of the fitting in a plot
                  k1k2="roy",                 # either "roy" (default, Roy1993a) or "lueker" (Lueker2000, calculated with seacarb) for K\_CO2 and K\_HCO3
                  khf="dickson",              # either "dickson" (default, Dickson1979a) or "perez" (Perez1987a, calculated with seacarb) for K\_HF}
                  datxbegin=0,                # at what x value (amount of titrant added) does the supplied curve start? (i.e. is the complete curve supplied or just a part?)
                  SumCO2Zero=FALSE)           # should SumCO2==0 ?
  {
    residuals <- function(state)
      {     
        if (K_CO2fit) {if (Evals){k_co2 <- state[[4]]/1e4} else {k_co2 <- state[[3]]/1e4}} # devide by 1e4 since the parameter was scaled to obtain parameters in the same range 

        if (SumCO2Zero) {state[[1]] <- 0} # should SumCO2==0?
        
        tit <- titration(aquaenv(S=ae$S, t=ae$t, p=ae$p,
                                 SumCO2=state[[1]], SumNH4=ae$SumNH4, SumH2S=ae$SumH2S, SumH3PO4=ae$SumH3PO4, SumSiOH4=ae$SumSiOH4, SumHNO3=ae$SumHNO3, SumHNO2=ae$SumHNO2,
                                 SumBOH3=ae$SumBOH3, SumH2SO4=ae$SumH2SO4, SumHF=ae$SumHF,
                                 TA=state[[2]], speciation=FALSE, skeleton=TRUE,
                                 k_w=k_w, k_co2=k_co2, k_hco3=k_hco3, k_boh3=k_boh3, k_hso4=k_hso4, k_hf=k_hf, k1k2=k1k2, khf=khf),
                         mass_sample=mass_sample, mass_titrant=titcurve[,1][[length(titcurve[,1])]], conc_titrant=conc_titrant,
                         S_titrant=S_titrant, steps=steps, type=type, seawater_titrant=seawater_titrant,
                         k_w=k_w, k_co2=k_co2, k_hco3=k_hco3, k_boh3=k_boh3, k_hso4=k_hso4, k_hf=k_hf, k1k2=k1k2, khf=khf) 
              
        calc <- tit$pH
       
        if (Evals)
          {
            # The Nernst equation relates E to TOTAL [H+] (DOE1994, p.7, ch.4, sop.3, Dickson2007)
            calc <- convert(calc, "pHscale", "free2tot", S=tit$S, t=tit$t, p=tit$p, SumH2SO4=tit$SumH2SO4, SumHF=tit$SumHF, khf=khf)
           
            # electrode polarity: E = E0 -(RT/F)ln([H+]) ("pos") or -E = E0 -(RT/F)ln([H+]) ("neg")
            if (electrode_polarity=="pos")
              {
                pol <- 1
              }
            else
              {
                pol <- -1
              }
            
            #Nernst Equation: E = E0 -(RT/F)ln([H+]); 83.14510 (bar*cm3)/(mol*K) = 8.314510 J/(mol*K): division by 10
            calc <- pol*(state[[3]]*1e2 - (((PhysChemConst$R/10)*ae$T)/PhysChemConst$F)*log(10^(-calc)))
           
            ylab <- "E/V"         
          }
        else if (!(pHscale=="free"))
          {
            calc <- convert(calc, "pHscale", paste("free2", pHscale, sep=""), S=tit$S, t=tit$t, p=tit$p, SumH2SO4=tit$SumH2SO4, SumHF=tit$SumHF, khf=khf) #conversion done here not on data before call, since S changes along the titration!

            ylab <- paste("pH (", pHscale, " scale)", sep="")         
          }
        else
          {
             ylab <- paste("pH (", pHscale, " scale)", sep="")         
          }
        if (!equalspaced)
          {
            # to cope with not equally spaced delta values for the masses (i.e. per step of the titration, NOT the same amount of titrant has been added),
            # we fit a function through our calculated pH curve (NOT the data)
            calcfun   <- approxfun(tit$delta_mass_titrant[(stepsbeforecurve+1):(steps+1)], calc[(stepsbeforecurve+1):(steps+1)], rule=2)

            residuals <- titcurve[,2]-calcfun(titcurve[,1])
          }
        else
          {
            residuals <- (titcurve[,2]-calc[(stepsbeforecurve+1):(steps+1)])
          }
        
        if (debug) # debug mode: make some variables global
          {
            tit  <<- tit
            calc <<- calc
          }
        if (verbose) # show the traject of the fitting procedure in a plot
          {
            ylim=range(titcurve[,2], calc)
            xlim=range(tit$delta_mass_titrant, titcurve[,1])
            plot(titcurve[,1], titcurve[,2], xlim=xlim, ylim=ylim, type="l", xlab="delta mass titrant", ylab=ylab)
            par(new=TRUE)
            plot(tit$delta_mass_titrant, calc, xlim=xlim, ylim=ylim, type="l", col="red", xlab="", ylab="")
          }
        return(residuals)
      }

    # deal with just partially supplied curves
    stepsbeforecurve <- 0
    steps <- (length(titcurve[,1])-1)
    if (!(datxbegin==0))
      {
        stepsbeforecurve <- (titcurve[,1][[1]] / ((titcurve[,1][[length(titcurve[,1])]] - titcurve[,1][[1]])/(length(titcurve[,1])-1)))
        steps            <- stepsbeforecurve  + (length(titcurve[,1])-1)
      } 
    
    if (Evals && K_CO2fit)
      {
        out    <- nls.lm(fn=residuals, par=c(TASumCO2guess, TASumCO2guess, E0guess/1e2, ae$K_CO2*1e4), control=nlscontrol) #multiply K_CO2 with 1e4 to bring the variables to fit into the same order of magnitude
        result <- list(out$par[[2]], out$par[[1]], out$par[[3]]*1e2, out$par[[4]]/1e4, out$deviance)                       #same thing for E0 but the other way round and with 1e2
        
        attr(result[[3]], "unit")     <- "V"
        attr(result[[4]], "unit")     <- "mol/kg-soln"
        attr(result[[4]], "pH scale") <- "free"
        names(result)                 <- c("TA", "SumCO2", "E0", "K_CO2", "sumofsquares")    
      }
    else if (Evals)
      { 
        out    <- nls.lm(fn=residuals, par=c(TASumCO2guess, TASumCO2guess, E0guess/1e2), control=nlscontrol)
        result <- list(out$par[[2]], out$par[[1]], out$par[[3]]*1e2, out$deviance)
        
        attr(result[[3]], "unit") <- "V"
        names(result)             <- c("TA", "SumCO2", "E0", "sumofsquares")    
      }
    else if (K_CO2fit)
      {    
        out    <- nls.lm(fn=residuals, par=c(TASumCO2guess, TASumCO2guess, ae$K_CO2*1e4), control=nlscontrol) #multiply K_CO2 with 1e4 to bring the variables to fit into the same order of magnitude
        result <- list(out$par[[2]], out$par[[1]], out$par[[3]]/1e4, out$deviance)
        
        attr(result[[3]], "unit")     <- "mol/kg-soln"
        attr(result[[3]], "pH scale") <- "free"
        names(result)                 <- c("TA", "SumCO2", "K_CO2", "sumofsquares")    
      }
    else
      {
        out    <- nls.lm(fn=residuals, par=c(TASumCO2guess, TASumCO2guess), control=nlscontrol)
        result <- list(out$par[[2]], out$par[[1]], out$deviance)
        
        names(result) <- c("TA","SumCO2","sumofsquares")
      }

    if (SumCO2Zero) {result[[2]] <- 0} # should SumCO2==0?
    
    attr(result[[1]], "unit")     <- "mol/kg-soln"
    attr(result[[2]], "unit")     <- "mol/kg-soln"    

    if (debug) # debug mode: make some variables global
      {
        out  <<- out
      }
   
    return(result)  # a list of up to five values ([TA] in mol/kg-solution, [SumCO2] in mol/kg-solution, E0 in V, K_CO2 in mol/kg-solution and on free scale, sum of the squared residuals)
  }
