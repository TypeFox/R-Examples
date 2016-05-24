#' C4 Photosynthesis
#' 
#' @description An implementation of the A-Ci curve for C4 plants, based on von Caemmerer et al. (2000)
#' @author Rhys Whitley
#' @param Ci Intercellular CO2 concentration (ppm)
#' @param PPFD Photosynthetic photon flux density (mu mol m-2 s-1)
#' @param Tleaf Leaf temperature (C)
#' @param VPMAX25 The maximum rate of PEP carboxylation (mu mol m-2 s-1)
#' @param Vcmax Maximum rate of carboxylation (mu mol m-2 s-1) (at 25C)
#' @param JMAX25 Maximum electron transport rate (at 25C)
#' @param Vpr PEP regeneration (mu mol m-2 s-1)
#' @param alpha Fraction of PSII activity in the bundle sheath (-)
#' @param gbs Bundle sheath conductance (mol m-2 s-1)
#' @param O2 Mesophyll O2 concentration
#' @param x Partitioning factor for electron transport
#' @param THETA Shape parameter of the non-rectangular hyperbola
#' @param Q10 T-dependence parameter for Michaelis-Menten coefficients.
#' @param RD0 Respiration at base temperature (RTEMP) (mu mol m-2 s-1)
#' @param RTEMP Base leaf temperature for respiration (C)
#' @param TBELOW Below this T, respiration is zero.
#' @param DAYRESP Fraction respiration in the light vs. that measured in the dark
#' @param Q10F T-dependence parameter of respiration 
#' @param FRM Fraction of dark respiration that is mesophyll respiration (Rm)
#' @param \dots Further arguments (currently ignored).
#' @details Note that the temperature response parameters have been hardwired in this function, and are based on von Caemmerer (2000).
#' 
#' Note that it is not (yet) possible to fit this curve to observations of photosynthesis (see \code{\link{fitaci}} to fit the C3 model of photosynthesis).
#' @references Caemmerer, S.V., 2000. Biochemical Models of Leaf Photosynthesis. Csiro Publishing.
#' @examples
#' # Simulate a C4 A-Ci curve. 
#' aci <- AciC4(Ci=seq(5,600, length=101))
#' with(aci, plot(Ci, ALEAF, type='l', ylim=c(0,max(ALEAF))))
#' @export
AciC4 <- function(Ci,
	PPFD=1500, 
	Tleaf = 25,
	VPMAX25=120, 
	JMAX25=400,
  Vcmax=60, 
	Vpr=80,         
	alpha=0.0,		  
	gbs=3e-3,		 
	O2=210,			 
	x=0.4, 			 
	THETA=0.7,   
	Q10 = 2.3,
	RD0=1, 
  RTEMP=25, 
  TBELOW=0, 
  DAYRESP=1, 
  Q10F=2, 
	FRM=0.5	,...){

	TK <- Tleaf+273.15
	
	# Temperature effects on Vcmax, Vpmax and Jmax (Massad et al. 2007)
  # This function returns value between 0 and 1.
  Arrhenius <- function(TK, Ea, Hd, DS){
        exp( Ea*(TK-298)/(298*.Rgas()*TK) ) * 
        		(1+exp( (298*DS-Hd)/(298*.Rgas()) )) / 
        		(1+exp( (TK*DS-Hd)/(TK*.Rgas()) )) 
  }
	
	# Half the reciprocal for Rubisco specificity (NOT CO2 compensation point)
	low_gammastar <- 1.93e-4
	
	# Michaelis-Menten coefficients for CO2 (Kc, mu mol mol-1) and O (Ko, mmol mol-1) and combined (K)
  Kc <- 650*Q10^((Tleaf-25)/10)
  Kp <- 80*Q10^((Tleaf-25)/10)
  Ko <- 450*Q10^((Tleaf-25)/10)
  K <- Kc*(1+O2/Ko)
          
  # T effects according to Massad et al. (2007)
	Vcmax <- Vcmax*Arrhenius(TK, 67294, 144568, 472)
	Vpmax <- VPMAX25*Arrhenius(TK, 70373, 117910, 376)
	Jmax <- JMAX25*Arrhenius(TK, 77900, 191929, 627)
	
	# Day leaf respiration, umol m-2 s-1
  if (Tleaf > TBELOW) {
      Rd <- RD0 * exp(Q10F * (Tleaf - RTEMP)) * DAYRESP
  } else {
      Rd <- 0.0
  }
	Rm <-  FRM*Rd
	
	# PEP carboxylation rate
	Vp <- pmin(Ci*Vpmax/(Ci+Kp),Vpr)
	
	# Quadratic solution for enzyme limited C4 assimilation
	a.c <- 1 - (alpha*Kc)/(0.047*Ko)
	b.c <- -( (Vp-Rm+gbs*Ci) + (Vcmax-Rd) + gbs*K + 
			alpha*low_gammastar/0.047*( low_gammastar*Vcmax+Rd*Kc/Ko ) )
	c.c <- (Vcmax-Rd)*(Vp-Rm+gbs*Ci) - (Vcmax*gbs*low_gammastar*O2 + Rd*gbs*K)
		
	A.enzyme <- (-b.c - sqrt(b.c^2 - 4*a.c*c.c)) / (2*a.c)

	# Non-rectangular hyperbola describing light effect on electron transport rate (J)
  Qp2 <- PPFD*(1-0.15)/2
  J <- (1/(2*THETA))*(Qp2+Jmax - sqrt((Qp2+Jmax)^2-4*THETA*Qp2*Jmax))
	
	# Quadratic solution for light-limited C4 assimilation
	a.j <- 1 - 7*low_gammastar*alpha/(3*0.047)
	b.j <- -( (x*J/2-Rm+gbs*Ci) + ((1-x)*J/3-Rd) + gbs*(7*low_gammastar*O2/3)
				+ alpha*low_gammastar/0.047*((1-x)*J/3+Rd) )
	c.j <- ( (x*J/2-Rm+gbs*Ci)*((1-x)*J/3-Rd) - gbs*low_gammastar*O2*((1-x)*J/3-7*Rd/3) )
		
	A.light <- (-b.j - sqrt(b.j^2 - 4*a.j*c.j)) / (2*a.j)

	# Actual assimilation rate
	An <- pmin(A.enzyme,A.light)
	Ac <- A.enzyme
	Aj <- A.light
	
	# Hyperbolic minimum (Buckley), to avoid discontinuity at transition from Ac to Aj
	shape2 <- 0.999
	Ad <- (Ac+Aj - sqrt((Ac+Aj)^2-4*shape2*Ac*Aj))/(2*shape2) - Rd
	Ac <- Ac - Rd
	Aj <- Aj - Rd
		
	return(data.frame(Ci=Ci, ALEAF=Ad, An=An, Ac=Ac, Aj=Aj, Vp=Vp, Rd=Rd, Tleaf=Tleaf, PPFD=PPFD))
}
