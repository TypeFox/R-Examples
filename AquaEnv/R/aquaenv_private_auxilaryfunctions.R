##################################################
# operations on objects of type AQUAENV
###################################################


# PRIVATE function:
# converts all elements of a special unit or pH scale in an object of class aquaenv                   
convert.aquaenv <- function(aquaenv,           # object of class aquaenv 
                            from,              # the unit which needs to be converted (as a string; must be a perfect match)
                            to,                # the unit to which the conversion should go
                            factor,            # the conversion factor to be applied: can either be a number (e.g. 1000 to convert from mol to mmol), or any of the conversion factors given in an object of class  aquaenv 
                            convattr="unit",   # which attribute should be converted? can either be "unit" or "pH scale"
                            ...) 
  {
    for (x in names(aquaenv))
      {
        if (!is.null(attr(aquaenv[[x]], convattr)))
            {
              if (attr(aquaenv[[x]], convattr) == from)
                {
                  aquaenv[[x]]                 <- aquaenv[[x]] * factor
                  attr(aquaenv[[x]], convattr) <- to
                }
            }
      }
    return(aquaenv)                            # object of class aquaenv whith the converted elements
  }


# PRIVATE function:
# returns the (maximal) length of the elements in an object of class aquaenv (i.e. > 1 if one of the input variables was a vector)
length.aquaenv <- function(x,           # object of class aquaenv
                           ...)
  {
    for (e in x)
      {
        if (length(e) > 1)
          {
            return(length(e))                 # the maximal length of the elements in the object of class aquaenv
          }
      }
    return(1)                                 # the maximal length of the elements in the object of class aquaenv
  }


# PRIVATE function:
# adds an element to an object of class aquaenv
c.aquaenv <- function(aquaenv,                # object of class aquaenv
                      x,                      # a vector of the form c(value, name) representing the element to be inserted into the object of class aquaenv
                      ...)
  {
    aquaenv[[x[[2]]]] <- x[[1]]
    return(aquaenv)                           # object of class aquaenv with the added element
  }


# PRIVATE function:
# merges the elements of two objects of class aquaenv: element names are taken from the first argument, the elements of which are also first in the merged object
merge.aquaenv <- function(x,           # object of class aquaenv: this is where the element names are taken from
                          y,           # object of class aquaenv: must contain at leas all the element (names) as x, extra elements are ignored
                          ...)
  {
    nam <- names(x)
    for (n in nam)
      {
        unit <- attr(x[[n]], "unit")
        x[[n]] <- c(x[[n]], y[[n]])
        attr(x[[n]], "unit") <- unit
      }
    return(x)                          # object of class aquaenv with merged elements
  }


# PRIVATE function:
# clones an object of class aquaenv: it is possible to supply a new value for either TA or pH; the switches speciation, skeleton, revelle, and dsa are obtained from the object to be cloned
cloneaquaenv <- function(aquaenv,             # object of class aquaenv
                         TA=NULL,             # optional new value for TA
                         pH=NULL,             # optional new value for pH
                         k_co2=NULL,          # used for TA fitting: give a K_CO2 and NOT calculate it from T and S: i.e. K_CO2 can be fitted in the routine as well
                         k1k2="roy",          # either "roy" (default, Roy1993a) or "lueker" (Lueker2000, calculated with seacarb) for K\_CO2 and K\_HCO3
                         khf="dickson",       # either "dickson" (default, Dickson1979a) or "perez" (Perez1987a, calculated with seacarb) for K\_HF}
                         khso4="dickson")     # either 'dickson" (default, Dickson1990) or "khoo" (Khoo1977) for K\_HSO4
  {
    if (is.null(TA) && is.null(pH))
      {
        pH <- aquaenv$pH
      }
    res <- aquaenv(S=aquaenv$S, t=aquaenv$t, p=aquaenv$p, SumCO2=aquaenv$SumCO2, SumNH4=aquaenv$SumNH4, SumH2S=aquaenv$SumH2S,SumH3PO4=aquaenv$SumH3PO4,
                   SumSiOH4=aquaenv$SumSiOH4, SumHNO3=aquaenv$SumHNO3, SumHNO2=aquaenv$SumHNO2, SumBOH3=aquaenv$SumBOH3, SumH2SO4=aquaenv$SumH2SO4,
                   SumHF=aquaenv$SumHF, pH=pH, TA=TA,
                   speciation=(!is.null(aquaenv$HCO3)), skeleton=(is.null(aquaenv$Na)), revelle=(!is.null(aquaenv$revelle)), dsa=(!is.null(aquaenv$dTAdH)), k1k2=k1k2, khf=khf, khso4=khso4)
    if (!is.null(k_co2))
      {
        res$K_CO2                   <- rep(k_co2,length(res))
        attr(res$K_CO2, "unit")     <- "mol/kg-soln"
        attr(res$K_CO2, "pH scale") <- "free"
      }
    
    return(res)                                   # cloned object of class aquaenv
  }


# PRIVATE function:
# creates an object of class aquaenv from a data frame (e.g. as supplied from the numerical solver of a dynamic model)
from.data.frame <- function(df)               # data frame
  {
    temp        <- as.list(df)
    class(temp) <- "aquaenv"
    return(temp)                              # object of class aquaenv
  }





#########################################################
# CONVERSION functions
#########################################################


# PRIVATE function:
# converts either the pH scale of a pH value, the pH scale of a dissociation constant (K*), or the unit of a concentration value
convert.standard <- function(x,               # the object to be converted (pH value, K* value, or concentration value)
                             vartype,         # the type of x, either "pHscale", "KHscale", or "conc"
                             what,            # the type of conversion to be done, for pH scales one of "free2tot", "free2sws", "free2nbs", ... (any combination of "free", "tot", "sws", and "nbs"); for concentrations one of "molar2molal", "molar2molin", ... (any combination of "molar" (mol/l), "molal" (mol/kg-H2O), and "molin" (mol/kg-solution))
                             S,               # salinity (in practical salinity units: no unit)
                             t,               # temperature in degrees centigrade 
                             p=0,             # gauge pressure (total pressure minus atmospheric pressure) in bars
                             SumH2SO4=NULL,   # total sulfate concentration in mol/kg-solution; if not supplied this is calculated from S
                             SumHF=NULL,      # total fluoride concentration in mol/kg-solution; if not supplied this is calculated from S
                             khf="dickson",   # either "dickson" (default, Dickson1979a) or "perez" (Perez1987a) for K\_HF
                             khso4="dickson") # either 'dickson" (default, Dickson1990) or "khoo" (Khoo1977) for K\_HSO4
  {
    result <- (switch
               (vartype,
                pHscale = switch
                (what,
                 free2tot = x - log10(scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$free2tot),
                 free2sws = x - log10(scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$free2sws),
                 free2nbs = x - log10(scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$free2nbs),
                 tot2free = x - log10(scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$tot2free),
                 tot2sws  = x - log10(scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$tot2sws),
                 tot2nbs  = x - log10(scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$tot2nbs),
                 sws2free = x - log10(scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$sws2free),
                 sws2tot  = x - log10(scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$sws2tot),
                 sws2nbs  = x - log10(scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$sws2nbs),
                 nbs2free = x - log10(scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$nbs2free),
                 nbs2tot  = x - log10(scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$nbs2tot),
                 nbs2sws  = x - log10(scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$nbs2sws)
                 ),
                KHscale = switch
                (what,
                 free2tot = x * scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$free2tot,
                 free2sws = x * scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$free2sws,
                 free2nbs = x * scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$free2nbs,
                 tot2free = x * scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$tot2free,
                 tot2sws  = x * scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$tot2sws,
                 tot2nbs  = x * scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$tot2nbs,
                 sws2free = x * scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$sws2free,
                 sws2tot  = x * scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$sws2tot,
                 sws2nbs  = x * scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$sws2nbs,
                 nbs2free = x * scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$nbs2free,
                 nbs2tot  = x * scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$nbs2tot,
                 nbs2sws  = x * scaleconvert(S, t, p, SumH2SO4, SumHF, khf, khso4)$nbs2sws
                 ),
                conc = switch
                (what,
                 molar2molal = x * (1/((seadensity(S,t)/1e3)* molal2molin(S))),
                 molar2molin = x * (1/((seadensity(S,t))/1e3))                ,
                 molal2molar = x * (molal2molin(S) * (seadensity(S,t))/1e3)   ,
                 molal2molin = x * (molal2molin(S))                            ,
                 molin2molar = x * (seadensity(S,t)/1e3)                      ,
                 molin2molal = x * (1/molal2molin(S))
                 )
                )
               )
    if ((what == "tot2free") || (what == "sws2free") || (what == "nbs2free"))
      {
        attr(result, "pH scale") <- "free"
      }
    else if ((what == "free2tot") || (what == "sws2tot") || (what == "nbs2tot"))
      {
        attr(result, "pH scale") <- "tot"
      }
    else if ((what == "free2nbs") || (what == "sws2nbs") || (what == "tot2nbs"))
      {
        attr(result, "pH scale") <- "nbs"
      }
    else if ((what == "molar2molal") || (what == "molin2molal"))
      {
        attr(result, "unit") <- "mol/kg-H2O"
      }
    else if ((what == "molal2molin") || (what == "molar2molin"))
      {
        attr(result, "unit") <- "mol/kg-soln"
      }
    else if ((what == "molal2molar") || (what == "molin2molar"))
      {
        attr(result, "unit") <- "mol/l-soln"
      }
    return(result)                            # converted pH, K*, or concentration value, attributed with the new unit/pH scale
  }


# PUBLIC function:
# calculates the depth (in m) from the gauge pressure p (or the total pressure P) and the latitude (in degrees: -90 to 90) and the atmospheric pressure (in bar)
# references Fofonoff1983
watdepth <- function(P=Pa, p=pmax(0, P-Pa), lat=0, Pa=1.013253) 
{
  gravity <- function(lat) 
    {
      X <- sin(lat * pi/180)
      X <- X * X
      grav = 9.780318 * (1 + (0.0052788 + 2.36e-05 * X) * X)
      return(grav)
    }

    P <- p * 10
    denom = gravity(lat) + 1.092e-06 * P
    nom = (9.72659 + (-2.2512e-05 + (2.279e-10 - 1.82e-15 * P) * 
        P) * P) * P
    return(nom/denom)
}


# PUBLIC function:
# calculates the gauge pressure from the depth (in m) and the latitude (in degrees: -90 to 90) and the atmospheric pressure (in bar) 
# references Fofonoff1983
gauge_p <- function(d, lat=0, Pa=1.01325)
  {
    gauge_p <- c()
    for (de in d)
      {
        if (de==0)
          {
            gauge_p <- c(gauge_p,0)
          }
        else
          {
            xx <- function(x)
              {
                return(de - watdepth(P=x+Pa, lat=lat))
              }
            gauge_p <- c(gauge_p, (uniroot(f=xx, interval=c(0,1300), tol=Technicals$uniroottol, maxiter=Technicals$maxiter)$root))
          }
      }
    return(gauge_p)
  }


# PRIVATE function:
# calculates the ionic strength I as a function of salinity S
# references: DOE1994, Zeebe2001, Roy1993b (the carbonic acid paper)
I <- function(S)                              # salinity S in practical salinity units (i.e. no unit)
  {
    return(19.924*S/(1000-1.005*S))           # ionic strength in mol/kg-solution (molinity)
  }


# PRIVATE function:
# calculates chlorinity Cl from salinity S
# references: DOE1994, Zeebe2001
Cl <- function(S)                             # salinity S in practical salinity units (i.e. no unit)
  {
    return(S/1.80655)                         # chlorinity Cl in permil
  }


# PRIVATE function:
# calculates concentrations of constituents of natural seawater from a given salinity S
# reference: DOE1994
seaconc <- function(spec,                     # constituent of seawater (chemical species) of which the concentration should be calculated. can be any name of the vectors ConcRelCl and MeanMolecularMass: "Cl", "SO4", "Br", "F", "Na", "Mg", "Ca", "K", "Sr", "B", "S"
                    S)                        # salinity S in practical salinity units (i.e. no unit)
  {
    return(                                   # concentration of the constituent of seawater speciefied in spec in mol/kg-solution (molinity): this is determined by the data in ConcRelCl and MeanMolecularMass
           ConcRelCl[[spec]]/MeanMolecularMass[[spec]]*Cl(S))   
  }
                    

# PRIVATE function:
# calculates the conversion factor converting from molality (mol/kg-H2O) to molinity (mol/kg-solution) from salinity S
# reference: Roy1993b (the carbonic acid paper), DOE1994
molal2molin <- function(S)                    # salinity S in practical salinity units (i.e. no unit)  
  {
    return(1-0.001005*S)                      # the conversion factor from molality (mol/kg-H2O) to molinity (mol/kg-solution)
  }


# PRIVATE function:
# calculates the temperature in Kelvin from the temperature in degrees centigrade
T <- function(t)                            # temperature in degrees centigrade
  {
    return(t - PhysChemConst$absZero)            # temperature in Kelvin
  }


# PRIVATE function:
# provides pH scale conversion factors (caution: the activity coefficient for H+ (needed for NBS scale conversions) is calculated with the Davies equation (Zeebe2001) which is only accurate up to ionic strengthes of I = 0.5)
# references: Dickson1984, DOE1994, Zeebe2001
scaleconvert <- function(S,                    # salinity S in practical salinity units (i.e. no unit)  
                         t,                    # temperature in degrees centigrade
                         p=0,                  # gauge pressure (total pressure minus atmospheric pressure) in bars
                         SumH2SO4=NULL,        # total sulfate concentration in mol/kg-solution; if not supplied this is calculated from S
                         SumHF=NULL,           # total fluoride concentration in mol/kg-solution; if not supplied this is calculated from S
                         khf="dickson",        # either "dickson" (Dickson1979a) or "perez" (Perez1987a) for K_HF
                         khso4="dickson")      # either 'dickson" (default, Dickson1990) or "khoo" (Khoo1977) for K\_HSO4
  {
    if (is.null(SumH2SO4))
      {
        SumH2SO4 = seaconc("SO4", S)
      }
    if (is.null(SumHF))
      {
        SumHF = seaconc("F", S)
      }

    K_HSO4 <- K_HSO4(S, t, p, khso4=khso4)
    K_HF   <- K_HF(S, t, p, SumH2SO4, SumHF, khf=khf)
    
    FreeToTot <- (1 + (SumH2SO4/K_HSO4))
    FreeToSWS <- (1 + (SumH2SO4/K_HSO4) + (SumHF/K_HF))
    attributes(FreeToTot) <- NULL
    attributes(FreeToSWS) <- NULL

    #davies equation: only valid up to I=0.5     
    SQRTI   <- sqrt(I(S))
    eT      <- PhysChemConst$e*T(t)
    A       <- 1.82e6/(eT*sqrt(eT))
    gamma_H <- 10^-((A*((SQRTI/(1+SQRTI)) - 0.2*I(S))))
      
    NBSToFree <- 1/(gamma_H*FreeToSWS) * molal2molin(S)      #Lewis1998, Perez1984: pH_NBS = -log10(gamma_H (H + HSO4 + HF))  with concs being molal
                                                             #i.e.: the NBS scale is related to the SEAWATER scale via gamma_H not the free scale
                                                             #      (since, if you measure with NBS buffers in seawater, you do not get the activity of the proton alone
                                                             #       but of the proton plus HSO4 and HF)
                                                             ##################
                                                             #Lewis1998: NBS scale is based on mol/kg-H2O (molality) and all other scales (incl free) on mol/kg-soln (molinity)
               
    return(list(                              # list of conversion factors "free2tot", "free2sws", etc.
                free2tot = FreeToTot,
                free2sws = FreeToSWS,
                free2nbs = 1/NBSToFree,
                tot2free = 1/FreeToTot,
                tot2sws  = 1/FreeToTot * FreeToSWS,
                tot2nbs  = 1/FreeToTot * 1/NBSToFree,
                sws2free = 1/FreeToSWS,
                sws2tot  = 1/FreeToSWS * FreeToTot,
                sws2nbs  = 1/FreeToSWS * 1/NBSToFree,
                nbs2free = NBSToFree,
                nbs2tot  = NBSToFree * FreeToTot,
                nbs2sws  = NBSToFree * FreeToSWS
                ))
  }


# PRIVATE function:
# calculates seawater density (in kg/m3) from temperature (in degrees centigrade) and salinity
# references: Millero1981, DOE1994
seadensity <- function(S,                    # salinity S in practical salinity units (i.e. no unit)  
                       t)                    # temperature in degrees centigrade

  {
    t2 <- t^2
    t3 <- t2 * t
    t4 <- t3 * t
    t5 <- t4 * t
            
    A <- 8.24493e-1    - 4.0899e-3*t + 7.6438e-5*t2 - 8.2467e-7*t3 + 5.3875e-9*t4
    B <- -5.72466e-3 + 1.0227e-4*t - 1.6546e-6*t2
    C <- 4.8314e-4
    
    densityWater <- 999.842594 + 6.793952e-2*t - 9.095290e-3*t2 + 1.001685e-4*t3 - 1.120083e-6*t4 + 6.536332e-9*t5
        
    return(densityWater + A*S + B*S*sqrt(S) + C*S^2)  # seawater density in kg/m        
  }





################################################################
# input / output (IO) functions
################################################################


# PRIVATE function:
# basic wrapper for the R plot function for plotting objects of class aquaenv; no return value, just side-effect
basicplot <- function(aquaenv,                # object of class aquaenv
                      xval,                   # x-value: the independent variable describing a change in elements of an object of class aquaenv
                      type="l",               # standard plot parameter;     default: plot lines
                      mgp=c(1.8, 0.5, 0),     # standard plot parameter;     default: axis title on line 1.8, axis labels on line 0.5, axis on line 0
                      mar=c(3,3,0.5,0.5),     # standard plot parameter;     default: margin of 3 lines bottom and left and 0.5 lines top and right
                      oma=c(0,0,0,0),         # standard plot parameter;     default: no outer margin
                      size=c(15,13),          # the size of the plot device; default: 15 (width) by 13 (height) inches
                      mfrow=c(11,10),         # standard plot parameter;     default: 11 columns and 10 rows of plots
                      device="x11",           # the device to plot on;       default: "x11" (can also be "eps" or "pdf")
                      filename="aquaenv",     # filename to be used if "eps" or "pdf" is selected for device
                      newdevice,              # flag: if TRUE, new plot device is opened
                      setpar,                 # flag: if TRUE parameters are set with the function par
		      ylab=NULL,              # y axis label: if given, it overrides the names from an aquaenv object
                      ...)
  {
    if (newdevice)
      {
        opendevice(device, size, filename)
      }
    if (setpar)
      {
        par(mfrow=mfrow, mar=mar, oma=oma, mgp=mgp)
      }
    aquaenv <- as.data.frame(aquaenv)
    for (i in 1:length(aquaenv))
      {
        if(is.null(ylab))
          {
            ylab_ <- names(aquaenv)[[i]]
          }
        else
          {
            ylab_ <- ylab
          }
        plot(xval, aquaenv[[i]], ylab=ylab_, type=type,  ...)
      }
  } 


# PRIVATE function:
# opens a device for plotting; no return value, just side-effect 
opendevice <- function(device,                 # either "x11", "eps", or "pdf"
                       size,                   # size of the plot device in the form c(width, height)
                       filename)               # filename to use if "eps" or "pdf" is used
  {
    if (device == "x11")
      {
        x11(width=size[[1]], height=size[[2]])
      }
    else if (device == "eps")
      {
        postscript(width=size[[1]], height=size[[2]], file=paste(filename, ".eps", sep=""), paper="special")
      }
    else if (device == "pdf")
      {
        pdf(width=size[[1]], height=size[[2]], file=paste(filename, ".pdf", sep=""), paper="special")
      }
  }


# PRIVATE function:
# plots all elements of an object of class aquaenv; no return value, just side-effect 
plotall <- function(aquaenv,                  # object of class aquaenv
                    xval,                     # x-value: the independent variable describing a change in elements of an object of class aquaenv
                    ...)
  {
    basicplot(aquaenv, xval=xval,  ...)
  }


# PRIVATE function:
# plots just the elements of an object of class aquaenv given in what; no return value, just side-effect 
selectplot <- function(aquaenv,               # object of class aquaenv
                       xval,                  # x-value: the independent variable describing a change in elements of an object of class aquaenv
                       what,                  # vector of names of elements of aquaenv that should be plotted
                       mfrow=c(1,1),          # standard plot parameter; default: just one plot
                       size=c(7,7),           # the size of the plot device; default: 7 (width) by 7 (height) inches
                       ...)
  {
    aquaenvnew <- aquaenv[what]
    class(aquaenvnew) <- "aquaenv"
    basicplot(aquaenvnew, xval=xval, mfrow=mfrow, size=size, ...)
  }


# PRIVATE function:
# creates a bjerrumplot from the elements of an object of class aquaenv given in what; no return value, just side-effect 
bjerrumplot <- function(aquaenv,              # object of class aquaenv
                        what,                 # vector of names of elements of aquaenv that should be plotted; if not specified:  what <- c("CO2", "HCO3", "CO3", "BOH3", "BOH4", "OH", "H3PO4", "H2PO4", "HPO4", "PO4", "SiOH4", "SiOOH3", "SiO2OH2", "H2S", "HS", "S2min", "NH4", "NH3", "H2SO4", "HSO4", "SO4", "HF", "F", "HNO3", "NO3", "HNO2", "NO2")
                        log=FALSE,            # should the plot be on a logarithmic y axis? 
                        palette=NULL,         # a vector of colors to use in the plot (either numbers or names given in colors())
                        device="x11",         # the device to plot on; default: "x11" (can also be "eps" or "pdf")
                        filename="aquaenv",   # filename to be used if "eps" or "pdf" is selected for device
                        size=c(12,10),        # the size of the plot device; default: 12 (width) by 10 (height) inches
                        ylim=NULL,            # standard plot parameter; if not supplied it will be calculated by range() of the elements to plot
                        lwd=2,                # standard plot parameter; width of the lines in the plot
                        xlab="free scale pH", # x axis label
                        mgp=c(1.8, 0.5, 0),   # standard plot parameter; default: axis title on line 1.8, axis labels on line 0.5, axis on line 0
                        mar=c(3,3,0.5,0.5),   # standard plot parameter; default: margin of 3 lines bottom and left and 0.5 lines top and right
                        oma=c(0,0,0,0),       # standard plot parameter; default: no outer margin
                        legendposition="bottomleft", # position of the legend
                        legendinset=0.05,     # standard legend parameter inset   
                        legendlwd=4,          # standard legend parameter lwd: line width of lines in legend
                        bg="white",           # standard legend parameter: default background color: white
                        newdevice,            # flag: if TRUE, new plot device is opened
                        setpar,               # flag: if TRUE parameters are set with the function par
                        ...)
  {
    if (is.null(what))
      {
        what <- c("CO2", "HCO3", "CO3", "BOH3", "BOH4", "OH", "H3PO4", "H2PO4", "HPO4", "PO4", "SiOH4", "SiOOH3", "SiO2OH2",
                   "H2S", "HS", "S2min", "NH4", "NH3", "H2SO4", "HSO4", "SO4", "HF", "F", "HNO3", "NO3", "HNO2", "NO2")
      }
    bjerrumvarslist <- aquaenv[what]
    class(bjerrumvarslist) <- "aquaenv"
    bjerrumvars <- as.data.frame(bjerrumvarslist)

    if (newdevice)
      {
        opendevice(device, size, filename)
      }
    if(setpar)
      {
        par(mar=mar, mgp=mgp, oma=oma)
      }
    
    if (is.null(palette))
      {
        palette <- 1:length(what)
      }

    if (log)
      {
        if (is.null(ylim))
          {
            ylim <- range(log10(bjerrumvars))
          }
        yvals <- log10(bjerrumvars)
        ylab  <- paste("log10([X]/(",attr(bjerrumvarslist[[1]], "unit"),"))", sep="")
      }
    else
      {
        if (is.null(ylim))
          {
            ylim  <- range(bjerrumvars)
          }
        yvals <- bjerrumvars
        ylab  <- attr(bjerrumvarslist[[1]], "unit")
      }
    for (i in 1:length(bjerrumvars))
      {
        plot(aquaenv$pH, yvals[[i]], type="l", ylab=ylab, xlab=xlab, ylim=ylim, col=palette[[i]], lwd=lwd, ...)
        par(new=TRUE)
      }
    par(new=FALSE)
    
    legend(legendposition, inset=legendinset, legend=names(bjerrumvarslist), col=palette, bg=bg, lwd=legendlwd, ...)
  }


# PRIVATE function:
# creates a cumulative plot from the elements of an object of class aquaenv given in what; no return value, just side-effect 
cumulativeplot <- function(aquaenv,           # object of class aquaenv
                           xval,              # x-value: the independent variable describing a change in elements of an object of class aquaenv
                           what,              # vector of names of elements of aquaenv that should be plotted
                           total=TRUE,        # should the sum of all elements specified in what be plotted as well?
                           palette=NULL,      # a vector of colors to use in the plot (either numbers or names given in colors())
                           device="x11",      # the device to plot on;       default: "x11" (can also be "eps" or "pdf")
                           filename="aquaenv",# filename to be used if "eps" or "pdf" is selected for device
                           size=c(12,10),     # the size of the plot device; default: 12 (width) by 10 (height) inches
                           ylim=NULL,         # standard plot parameter; if not supplied it will be calculated by an adaptation of range() of the elements to plot
                           lwd=2,             # standard plot parameter; width of the lines in the plot
                           mgp=c(1.8, 0.5, 0),# standard plot parameter; default: axis title on line 1.8, axis labels on line 0.5, axis on line 0
                           mar=c(3,3,0.5,0.5),# standard plot parameter; default: margin of 3 lines bottom and left and 0.5 lines top and right
                           oma=c(0,0,0,0),    # standard plot parameter; default: no outer margin
                           legendposition="bottomleft", # position of the legend
                           legendinset=0.05,  # standard legend parameter inset  
                           legendlwd=4,       # standard legend parameter lwd: line width of lines in legend
                           bg="white",        # standard legend parameter: default background color: white
                           y.intersp=1.2,     # standard legend parameter; default: 1.2 lines space between the lines in the legend
                           newdevice,         # flag: if TRUE, new plot device is opened
                           setpar,            # flag: if TRUE parameters are set with the function par
                           ...)
  {
    if (is.null(what))
      {
        what=names(aquaenv)
      }

    cumulativevarslist <- aquaenv[what]
    class(cumulativevarslist) <- "aquaenv"
    cumulativevars <- as.data.frame(cumulativevarslist)
    
    if (is.null(ylim))
      {
        ylim <- c(0,0)
        for (var in cumulativevars)
          {
            ylim <- ylim + range(c(var,0))
          }
      }
    if (is.null(palette))
      {
        palette <- 1:length(names(cumulativevars))
      }

    if (newdevice)
      {
        opendevice(device, size, filename)
      }
    if (setpar)
      {
        par(mar=mar, mgp=mgp, oma=oma)
      }
    plot(xval, rep(0,length(xval)), type="l", ylim=ylim, col="white", ...)
    sumfuncpos <- rep(0,length(xval))
    for (x in 1:(length(cumulativevars)))
      {
        yval <- (cumulativevars[[x]]) 
        yval[yval<=0] <- 0
        newsumfuncpos <- sumfuncpos + yval

        if (!(identical(yval, rep(0,length(yval)))))
          {
            polygon(c(xval,rev(xval)), c(newsumfuncpos, rev(sumfuncpos)), col=palette[[x]], border=NA)
          }
        
        sumfuncpos <- newsumfuncpos
      }
    sumfuncneg <- rep(0,length(xval))
    for (x in 1:(length(cumulativevars)))
      {
        yval <- (cumulativevars[[x]])
        yval[yval>=0] <- 0
        newsumfuncneg <- sumfuncneg + yval
        
        if (!(identical(yval, rep(0,length(yval)))))
          {
            polygon(c(xval,rev(xval)), c(newsumfuncneg, rev(sumfuncneg)), col=palette[[x]], border=NA)
          }
        
        sumfuncneg <- newsumfuncneg
      }
      if (total)
      {
        par(new=TRUE)
        plot(xval, apply(cumulativevars, 1, sum), col="gray", type="l", ylim=ylim, xlab="", ylab="", lwd=lwd)
      }	
    legend(legendposition, legend=names(cumulativevars), col=palette, inset=legendinset, y.intersp=y.intersp, bg=bg, lwd=legendlwd)    
  }

