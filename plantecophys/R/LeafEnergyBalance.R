#' Coupled leaf gas exchange model with energy balance
#' @description As \code{\link{Photosyn}}, but calculates the leaf temperature based on the leaf's energy balance. Including sensible and long-wave heat loss, latent heat loss from evaporation, and solar radiation input. 
#' 
#'@details Uses the Penman-Monteith equation to calculate the leaf transpiration rate, and finds Tleaf by solving the leaf energy balance iteratively. In the solution, it is accounted for that stomatal conductance (via the dependence of photosynthesis on Tleaf) and net radiation depend on Tleaf.
#'
#'Also included is the function \code{FindTleaf}, which calculates the leaf temperature if the stomatal conductance is known.
#'@param Tair Air temperature (C)
#'@param VPD The vapour pressure deficit of the air (i.e. not the leaf-to-air VPD) (kPa).
#'@param Wind Wind speed (m s-1)
#'@param Wleaf Leaf width (m)
#'@param StomatalRatio The stomatal ratio (cf. Licor6400 terminology), if it is 1, leaves have stomata only on one side (hypostomatous), 2 for leaves with stomata on both sides (amphistomatous).
#'@param LeafAbs Leaf absorptance of solar radiation (0-1).
#'@param RH The relative humidity of the air (i.e. not calculated with leaf temperature) (in percent).
#'@param gs For \code{FindTleaf}, the stomatal conductance (mol m-2 s-1).
#'@param \dots Further parameters passed to \code{\link{Photosyn}}. Note that Tleaf is not allowed as an input, since that is calculated by \code{PhotosynEB} from energy balance.
#'@export PhotosynEB
#'@rdname PhotosynEB
PhotosynEB <- function(Tair=25,
                       VPD=1.5,
                       Wind = 2,   
                       Wleaf = 0.02, 
                       StomatalRatio = 1,   # 2 for amphistomatous
                       LeafAbs = 0.86, 
                       RH=NULL,
                       ...){
  
  
  if(!is.null(RH))VPD <- RHtoVPD(RH,Tair)
  
  # Non-vectorized function declared here.
  PhotosynEBfun <- function(Tair,
                            VPD,
                            Wind,   
                            Wleaf, 
                            StomatalRatio,   # 2 for amphistomatous
                            LeafAbs,
                            ...){   # passed to Photosyn
    
    
    m <- match.call(expand.dots=TRUE)
    if("Tleaf" %in% names(m))
      stop("Cannot pass Tleaf to PhotosynEB - it is calculated from energy balance.")
    
    # Wrapper for Photosyn to find gs only  
    gsfun <- function(...)Photosyn(...)$GS
    
    # Find Tleaf. Here, we take into account that Tleaf as solved from
    # energy balance affects gs, so this is the second loop to solve for Tleaf.
    fx <- function(x, Tair, Wind, VPD, Wleaf, StomatalRatio, LeafAbs,  ...){
      newx <- FindTleaf(Tair=Tair, gs=gsfun(Tleaf=x, VPD=VPD, ...), 
                        Wind=Wind, Wleaf=Wleaf, 
                        StomatalRatio=StomatalRatio, LeafAbs=LeafAbs)
      (newx - x)^2
    }

    Tleaf <- try(optimize(fx, interval=c(max(1, Tair-15), Tair+15), Tair=Tair, Wind=Wind, Wleaf=Wleaf, 
                     VPD=VPD, StomatalRatio=StomatalRatio, LeafAbs=LeafAbs, ...)$minimum)
    
    failed <- FALSE
    if(inherits(Tleaf, "try-error")){
      Tleaf <- Tair
      failed <- TRUE
    }

    # Now run Photosyn
    p <- Photosyn(Tleaf=Tleaf, VPD=VPD, ...)
    
    # And energy balance components
    e <- LeafEnergyBalance(Tleaf=Tleaf, Tair=Tair, gs=p$GS, 
                           PPFD=p$PPFD, VPD=p$VPD, Patm=p$Patm, 
                           Wind=Wind, Wleaf=Wleaf, 
                           StomatalRatio=StomatalRatio, LeafAbs=LeafAbs, 
                           returnwhat="fluxes")
    res <- cbind(p,e)
    
    # Replace ELEAF with energy-balance one.
    res$ELEAF <- res$ELEAFeb
    res$ELEAFeb <- NULL
    res$VPDleaf <- VPDairToLeaf(res$VPD, Tair, res$Tleaf)
    res$Tair <- Tair
    res$Wind <- Wind
    res$failed <- failed
    
    return(res)
  }
  
  m <- t(mapply(PhotosynEBfun, Tair=Tair, 
                VPD=VPD, 
                Wind=Wind, Wleaf=Wleaf, StomatalRatio=StomatalRatio, 
                LeafAbs=LeafAbs,  ..., SIMPLIFY=FALSE))
  
  
  return(as.data.frame(do.call(rbind, m)))
}

# The net leaf energy balance, given that we know Tleaf, gs
LeafEnergyBalance <- function(Tleaf = 21.5, Tair = 20, 
                              gs = 0.15,
                              PPFD = 1500, VPD = 2, Patm = 101,
                              Wind = 2, Wleaf = 0.02, 
                              StomatalRatio = 1,   # 2 for amphistomatous
                              LeafAbs = 0.5,   # in shortwave range, much less than PAR
                              returnwhat=c("balance","fluxes")){   
                              
  
  returnwhat <- match.arg(returnwhat)
  
  # Constants
  Boltz <- 5.67 * 10^-8     # w M-2 K-4
  Emissivity <- 0.95        # -
  LatEvap <- 2.54           # MJ kg-1
  CPAIR <- 1010.0           # J kg-1 K-1
  
  H2OLV0 <- 2.501e6         # J kg-1
  H2OMW <- 18e-3            # J kg-1
  AIRMA <- 29.e-3           # mol mass air (kg/mol)
  AIRDENS <- 1.204          # kg m-3
  UMOLPERJ <- 4.57
  DHEAT <- 21.5e-6          # molecular diffusivity for heat

  
  
  # Density of dry air
  AIRDENS <- Patm*1000/(287.058 * Tk(Tair))
  
  # Latent heat of water vapour at air temperature (J mol-1)
  LHV <- (H2OLV0 - 2.365E3 * Tair) * H2OMW
  
  # Const s in Penman-Monteith equation  (Pa K-1)
  SLOPE <- (esat(Tair + 0.1) - esat(Tair)) / 0.1
  
  # Radiation conductance (mol m-2 s-1)
  Gradiation <- 4.*Boltz*Tk(Tair)^3 * Emissivity / (CPAIR * AIRMA)

  # See Leuning et al (1995) PC&E 18:1183-1200 Appendix E
  # Boundary layer conductance for heat - single sided, forced convection
  CMOLAR <- Patm*1000 / (8.314 * Tk(Tair))   # .Rgas() in package...
  Gbhforced <- 0.003 * sqrt(Wind/Wleaf) * CMOLAR
  
  # Free convection
  GRASHOF <- 1.6E8 * abs(Tleaf-Tair) * (Wleaf^3) # Grashof number
  Gbhfree <- 0.5 * DHEAT * (GRASHOF^0.25) / Wleaf * CMOLAR
  
  # Total conductance to heat (both leaf sides)
  Gbh <- 2*(Gbhfree + Gbhforced)
  
  # Heat and radiative conductance
  Gbhr <- Gbh + 2*Gradiation
  
  # Boundary layer conductance for water (mol m-2 s-1)
  Gbw <- StomatalRatio * 1.075 * Gbh  # Leuning 1995
  gw <- gs*Gbw/(gs + Gbw)
  
  # Longwave radiation
  # (positive flux is heat loss from leaf)
  Rlongup <- Emissivity*Boltz*Tk(Tleaf)^4
  
  # Rnet
  Rsol <- 2*PPFD/UMOLPERJ   # W m-2
  Rnet <- LeafAbs*Rsol - Rlongup   # full
  
  # Isothermal net radiation (Leuning et al. 1995, Appendix)
  ea <- esat(Tair) - 1000*VPD
  ema <- 0.642*(ea/Tk(Tair))^(1/7)
  Rnetiso <- LeafAbs*Rsol - (1 - ema)*Boltz*Tk(Tair)^4 # isothermal net radiation
  
  # Isothermal version of the Penmon-Monteith equation
  GAMMA <- CPAIR*AIRMA*Patm*1000/LHV
  ET <- (1/LHV) * (SLOPE * Rnetiso + 1000*VPD * Gbh * CPAIR * AIRMA) / (SLOPE + GAMMA * Gbhr/gw)

  # Latent heat loss
  lambdaET <- LHV * ET

  # Heat flux calculated using Gradiation (Leuning 1995, Eq. 11)
  Y <- 1/(1 + Gradiation/Gbh)
  H2 <- Y*(Rnetiso - lambdaET)
  
  # Heat flux calculated from leaf-air T difference.
  # (positive flux is heat loss from leaf)
  H <- -CPAIR * AIRDENS * (Gbh/CMOLAR) * (Tair - Tleaf)
  
  # Leaf-air temperature difference recalculated from energy balance.
  # (same equation as above!)
  Tleaf2 <- Tair + H2/(CPAIR * AIRDENS * (Gbh/CMOLAR))
  
  # Difference between input Tleaf and calculated, this will be minimized.
  EnergyBal <- Tleaf - Tleaf2
  
  if(returnwhat == "balance")return(EnergyBal)
  
  if(returnwhat == "fluxes"){
    
    l <- data.frame(ELEAFeb=1000*ET, Gradiation=Gradiation, Rsol=Rsol, Rnetiso=Rnetiso, Rlongup=Rlongup, H=H, lambdaET=lambdaET, gw=gw, Gbh=Gbh, H2=H2, Tleaf2=Tleaf2)
    return(l)
  }
}
  
  
#' Calculate Tleaf from energy balance, given that we know gs
#' @export
#' @rdname PhotosynEB
#' @importFrom stats uniroot
FindTleaf <- function(gs, Tair, ...){
  
 Tleaf <- try(uniroot(LeafEnergyBalance, interval=c(Tair-15, Tair+15), 
                    gs=gs, Tair=Tair, ...)$root)
#   Tleaf <- optimize(LeafEnergyBalance, interval=c(Tair-15, Tair+15), 
#                       gs=gs, Tair=Tair,...)$minimum
  
return(Tleaf)
}








