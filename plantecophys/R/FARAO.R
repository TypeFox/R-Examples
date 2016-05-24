#' FARquhar And Opti
#' @description The numerical solution of the optimal stomatal conductance model, coupled with the Farquhar model of photosynthesis. The model of Medlyn et al. (2011) is an approximation to this full numeric solution.
#' @param lambda The marginal cost of water (mol mol-1)
#' @param Ca The CO2 concentration. 
#' @param VPD Vapor pressure deficit (kPa)
#' @param photo Which photosynthesis rate should stomata respond to? Defaults to 'BOTH', i.e. the minimum of Vcmax and Jmax limited rates.
#' @param energybalance If TRUE (Default = FALSE), calculates leaf temperature from energy balance (and its effects on photosynthesis as well as leaf transpiration), using \code{\link{PhotosynEB}}.
#' @param C4 If TRUE, uses the C4 photosynthesis routine (\code{\link{AciC4}})
#' @param Tair Air temperature (deg C)
#' @param Wind Wind speed (m s-1) (only used if energybalance=TRUE)
#' @param Wleaf Leaf width (m) (only used if energybalance=TRUE)
#' @param StomatalRatio The stomatal ratio (see \code{\link{PhotosynEB}}) (only used if energybalance=TRUE)
#' @param LeafAbs Leaf absorptance (see \code{\link{PhotosynEB}}) (only used if energybalance=TRUE)
#' @param ... All other parameters are passed to \code{\link{Aci}}
#' @author Remko Duursma
#' @details This model finds the Ci that maximizes A - lambda*E (Cowan & Farquhar 1977, see also Medlyn et al. 2011). The new function FARAO2 is a much simpler (and probably more stable) implementation, based on Buckley et al. 2014 (P,C&E). Both functions are provided, as FARAO has a few more options than FARAO2, at the moment.
#' @references 
#' Buckley, T.N., Martorell, S., Diaz-Espejo, A., Tomas, M., Medrano, H., 2014. Is stomatal conductance optimized over both time and space in plant crowns? A field test in grapevine (Vitis vinifera). Plant Cell Environ doi:10.1111/pce.12343
#' 
#' Cowan, I. and G.D. Farquhar. 1977. Stomatal function in relation to leaf metabolism and environment. Symposia of the Society for Experimental Biology. 31:471-505.
#' 
#' Medlyn, B.E., R.A. Duursma, D. Eamus, D.S. Ellsworth, I.C. Prentice, C.V.M. Barton, K.Y. Crous, P. De Angelis, M. Freeman and L. Wingate. 2011. Reconciling the optimal and empirical approaches to modelling stomatal conductance. Global Change Biology. 17:2134-2144.
#' @export
#' @importFrom stats optimize
#' @rdname FARAO
FARAO <- function(lambda=0.002, Ca=400, VPD=1, 
                  photo=c("BOTH","VCMAX","JMAX"), 
                  energybalance=FALSE, 
                  C4=FALSE,                   
                  Tair=25,
                  Wind=2,
                  Wleaf=0.02,
                  StomatalRatio=1,
                  LeafAbs=0.86,
                  ...){
  
  photo <- match.arg(photo)
  
  # Non-vectorized form.
  FARAOfun <- function(lambda, Ca, VPD, 
                       photo, 
                       energybalance, C4,                   
                       Tair,
                       Wind,
                       Wleaf,
                       StomatalRatio,
                       LeafAbs,
                       ...){
    
  
    if(!energybalance){
      fx <- function(Ca,...)optimize(OPTfun, interval=c(0,Ca), 
                                     maximum=TRUE,Ca=Ca,
                                     ...)$maximum
      optimalcis <- mapply(fx,Ca=Ca,lambda=lambda, photo=photo, C4=C4, VPD=VPD,...)
      
      res <- as.data.frame(OPTfun(Ci=optimalcis, retobjfun=FALSE, 
                                Ca=Ca, photo=photo, C4=C4, VPD=VPD,...))
    } else {
      fx <- function(Ca,...)optimize(OPTfunEB, interval=c(0,Ca), 
                                     maximum=TRUE,Ca=Ca,
                                     ...)$maximum
      optimalcis <- mapply(fx,Ca=Ca,lambda=lambda, photo=photo, C4=C4, VPD=VPD,
                           Tair=Tair,
                           Wind=Wind,
                           Wleaf=Wleaf,
                           StomatalRatio=StomatalRatio,
                           LeafAbs=LeafAbs,
                           ...)
      
      res <- as.data.frame(OPTfunEB(Ci=optimalcis, retobjfun=FALSE, 
                                  Ca=Ca, photo=photo, C4=C4, VPD=VPD,
                                  Tair=Tair,
                                  Wind=Wind,
                                  Wleaf=Wleaf,
                                  StomatalRatio=StomatalRatio,
                                  LeafAbs=LeafAbs,...))
    }
    
    return(res)
  }
  
  f <- t(mapply(FARAOfun, lambda=lambda, Ca=Ca, VPD=VPD, 
              photo=photo, 
              energybalance=energybalance, C4=C4,
              Tair=Tair,
              Wind=Wind,
              Wleaf=Wleaf,
              StomatalRatio=StomatalRatio,
              LeafAbs=LeafAbs, ..., SIMPLIFY=FALSE))
  
return(as.data.frame(do.call(rbind,f)))

}

# This function returns the 'objective function' A - lambda*E
# This is to be optimized by the next function by varying ci.
OPTfun <- function(Ci,              # mu mol mol-1
                    lambda=0.002,   # mol mol-1
                    Ca=400,         # mu mol mol-1
                    VPD=1,        # kPa
					          Patm=101,         # ambient pressure, kPa
					          photo=c("BOTH","VCMAX","JMAX"),
                    energybalance=FALSE,
                    retobjfun=TRUE, # if false, returns A, g and E (otherwise sum(A-lambda*E))
					          C4=FALSE,
					          ...){     

  GCtoGW <- 1.57
	VPDmol <- VPD/Patm
	
	photo <- match.arg(photo)
	
  # Given a Ci, calculate photosynthetic rate
  if(!C4)
		run <- Aci(Ci=Ci, VPD=VPD, ...)   # note that VPD does not do anything, just for consistency in I/O
	else 
		run <- AciC4(Ci, VPD=VPD, ...)
	
  if(photo == "BOTH")A <- run$ALEAF
	if(photo == "VCMAX")A <- run$Ac
	if(photo == "JMAX")A <- run$Aj
  
  # Given Ci and A, calculate gs (diffusion constraint)
  gs <- GCtoGW * A / (Ca - Ci)
	    
  # Transpiration rate
  E <- gs*VPDmol
    
    
  # Objective function to be maximized (Cowan-Farquhar condition)
  objfun <- 10^-6*A - lambda*E

if(retobjfun)return(objfun)

if(!retobjfun)return(list( Ci=Ci, ALEAF=A, GS=gs, ELEAF=E*1000, Ac=run$Ac, Aj=run$Aj,
                           Rd=run$Rd, VPD=VPD, Tleaf=run$Tleaf,  Ca=Ca, PPFD=run$PPFD ))
}         


OPTfunEB <- function(Ci,           # mu mol mol-1
                   lambda=0.002,   # mol mol-1
                   Ca=400,         # mu mol mol-1
                   VPD=1.5,        # AIR VPD! kPa
                   Patm=101,       # ambient pressure, kPa
                   Tair=25,
                   Wind=2,
                   Wleaf=0.02,
                   StomatalRatio=1,
                   LeafAbs=0.86,
                   photo=c("BOTH","VCMAX","JMAX"),
                   retobjfun=TRUE, # if false, returns A, g and E (otherwise sum(A-lambda*E))
                   C4=FALSE,
                   ...){     
  
  GCtoGW <- 1.57
  VPDmol <- VPD/Patm
  
  photo <- match.arg(photo)
  
  gsfun <- function(Ci, VPD, returnwhat=c("gs","all"), ...){
    
    returnwhat <- match.arg(returnwhat)
    
    # Given a Ci, calculate photosynthetic rate
    if(!C4)
      run <- Aci(Ci, VPD=VPD, ...)   # note that VPD does not do anything, just for consistency in I/O
    else 
      run <- AciC4(Ci, VPD=VPD, ...)
    
    if(photo == "BOTH")A <- run$ALEAF
    if(photo == "VCMAX")A <- run$Ac
    if(photo == "JMAX")A <- run$Aj
    
    # Given Ci and A, calculate gs (diffusion constraint)
    gs <- GCtoGW * A / (Ca - Ci)
    
  if(returnwhat == "gs")return(gs)
  if(returnwhat == "all")return(list(run=run, GS=gs, A=A))
  }
  
  # Find Tleaf. Here, we take into account that Tleaf as solved from
  # energy balance affects gs (because it affects)
  fx <- function(x, Ci, Tair, Wind, VPD, Wleaf, StomatalRatio, LeafAbs, ...){
    newx <- FindTleaf(Tair=Tair, gs=gsfun(Ci=Ci, Tleaf=x, VPD=VPD, ...), 
                      Wind=Wind, Wleaf=Wleaf, 
                      StomatalRatio=StomatalRatio, LeafAbs=LeafAbs)
    (newx - x)^2
  }

  Tleaf <- optimize(fx, interval=c(Tair-10, Tair+10), Ci=Ci, Tair=Tair, Wind=Wind, VPD=VPD, Wleaf=Wleaf, 
                   StomatalRatio=StomatalRatio, LeafAbs=LeafAbs, ...)$minimum
  
  z <- gsfun(Ci=Ci, Tleaf=Tleaf, VPD=VPD, returnwhat="all",...)
  GS <- z$GS
  A <- z$A
  
  # And energy balance components
  e <- LeafEnergyBalance(Tleaf=Tleaf, Tair=Tair, gs=GS, 
                         PPFD=z$run$PPFD, VPD=VPD, Patm=z$run$Patm, 
                         Wind=Wind, Wleaf=Wleaf, 
                         StomatalRatio=StomatalRatio, LeafAbs=LeafAbs,
                         returnwhat="fluxes")

  E <- e$ELEAFeb
  
  # Objective function to be maximized (Cowan-Farquhar)
  objfun <- 10^-6*A - lambda*E/1000
  
  if(retobjfun)return(objfun)
  
  if(!retobjfun)return(list( Ci=Ci, ALEAF=A, GS=GS, ELEAF=E, Ac=z$run$Ac, Aj=z$run$Aj,
                             Rd=z$run$Rd, VPD=VPD, Tleaf=Tleaf,  Ca=Ca, PPFD=z$run$PPFD ))
}  


getdAdE <- function(Ci,...,energybalance=FALSE,
                       returnwhat=c("dAdE","both")){
  
  returnwhat <- match.arg(returnwhat)
  delta <- 1e-03
  
  if(energybalance){
    r1 <- PhotosynEB(Ci=Ci, ...)
    r2 <- PhotosynEB(Ci=Ci+delta, ...)
  } else {
    r1 <- Photosyn(Ci=Ci, ...)
    r2 <- Photosyn(Ci=Ci+delta, ...)
  }
    
  dA <- r2$ALEAF - r1$ALEAF
  dE <- r2$ELEAF - r1$ELEAF
  
  if(returnwhat == "dAdE")
    return(dA/dE)
  else
    return(c(dA=dA, dE=dE))
  
}

#' @rdname FARAO
#' @export
FARAO2 <- function(lambda=0.002, Ca=400, energybalance=FALSE, ...){
  
  
  faraofun <- function(lambda,Ca,energybalance,...){
    f <- function(x, ...)(getdAdE(x, energybalance=energybalance, ...) - lambda*1000)^2
    
    CI <- try(optimize(f, c(80, Ca-0.1), ...)$minimum)
    if(inherits(CI, "try-error"))return(NULL)
    
    if(energybalance)
      p <- PhotosynEB(Ci=CI, ...)
    else
      p <- Photosyn(Ci=CI, ...)
  
  return(p)
  }
  
  m <- mapply(faraofun, lambda=lambda, Ca=Ca, energybalance=energybalance, ..., SIMPLIFY=FALSE)
  
  
return(do.call(rbind,m))
}


