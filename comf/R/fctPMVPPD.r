# Functions return:
#
# pmv - Predicted mean Vote
# ppd - Predicted Percentage of Dissatisfied
#
# File contains 1 function:
#   - calcPMVPPD(ta, tr, vel, rh, clo=.5, met=1, wme=0)
#       returns pmv, ppd
#
# v1.0 done by Marcel Schweiker 


# Function: pmv/ppd by Fanger ########################
###########################################
calcPMVPPD <- function(ta, tr, vel, rh, clo=.5, met=1, wme=0, basMet=58.15){
  
  m   <- met * basMet 
  w   <- wme * basMet
  mw  <- m - w
  icl <- .155 * clo
  pa  <- rh * 10 * exp(16.6536 - (4030.183 / (ta + 235)))
  
  # Compute the corresponding fcl value
  if (icl <= .078){
    fcl <- 1 + 1.29 * icl
  } else {
    fcl <- 1.05 + .645 * icl
  }
  
  fcic <- icl * fcl
  p2   <- fcic * 3.96
  p3   <- fcic * 100
  tra  <- tr + 273
  taa  <- ta + 273
  p1   <- fcic * taa
  p4   <- 308.7 - .028 * mw + p2 * (tra / 100) ^ 4
  
  # First guess for surface temperature
  tclA <- taa + (35.5-ta) / (3.5 * (6.45 * icl + .1))
  xn   <- tclA / 100
  xf   <- xn
  hcf  <- 12.1 * (vel) ^ .5
  noi  <- 0
  eps  <- .00015
  
  
  # COmPUTE sURFAce TEmPEratuRE OF cloTHING BY ITEraTIONs
  while (noi < 150){
    xf  <- (xf + xn) / 2
    hcn <- 2.38 * abs(100 * xf - taa) ^ .25
    if (hcf > hcn){
      hc <- hcf
    } else {
      hc <- hcn
    }
    xn  <- (p4 + p1 * hc - p2 * xf ^ 4) / (100 + p3 * hc)
    noi <- noi + 1
    if(noi > 1 & abs(xn - xf) <= eps){break} 
  }
  tcl <- 100 * xn - 273
  
  # COmPUTE pmv

  pm1 <- 3.96 * fcl * (xn ^ 4 - (tra / 100) ^ 4)
  pm2 <- fcl * hc * (tcl - ta)
  pm3 <- .303 * exp(-.036 * m) + .028
  if (mw > basMet){
     pm4 <- .42 * (mw - basMet) 
  } else {
    pm4 <- 0
  }
  pm5 <- 3.05 * .001 * (5733 - 6.99 * mw - pa)
  pm6 <- 1.7 * .00001 * m * (5867 - pa) + .0014 * m * (34 - ta)
  pmv <- pm3 * (mw - pm5 - pm4 - pm6 - pm1 - pm2)
    
  ppd <- 100 - 95 * exp(-.03353 * pmv ^ 4 - .2179 * pmv ^ 2)
  data.frame(pmv, ppd)
}
# End calcpmv ##############################
###########################################

## calc apmv ################################
## based on Yao et al.
############################################
calcaPMV <- function(ta, tr, vel, rh, clo = .5, met = 1, wme = 0, apCoeff){

	pmv  <- calcPMVPPD(ta, tr, vel, rh, clo, met, wme)[1]
	apmv <- pmv / (1 + apCoeff * pmv)
	names(apmv) <- "apmv"
	data.frame(apmv = apmv)

}

## calc epmv ################################
## based on Fanger et al.
############################################
calcePMV <- function(ta, tr, vel, rh, clo = .5, met = 1, wme = 0, epCoeff, asv){

	met <- ifelse (asv > 0, met * (1 + asv * (-.067)), met)
	pmv <- calcPMVPPD(ta, tr, vel, rh, clo, met, wme)[1]
	epmv <- epCoeff * pmv
	names(epmv) <- "epmv"
	data.frame(epmv = epmv)

}


## calc pmv adj ##############################
## based on ASHRAE 55-2013
#############################################

calcpmvadj <- function(ta, tr, vel, rh, clo, met, wme = 0){

f <- function(x){calc2Node(ta, tr, vel, rh, clo, met, wme, pb = 760, ltime = 60, ht = 171, wt = 70, tu = 40, obj = "pmvadj")[2] - calc2Node(ta + x, tr + x, .1, rh, clo, met, wme, pb = 760, ltime = 60, ht = 171, wt = 70, tu = 40, obj = "pmvadj")[2]}
ce <- bisect(f, -15, 1)$x

pmvadj <- calcPMVPPD(ta + ce, tr + ce, .1, rh, clo, met, wme)[1]
names(pmvadj) <- "pmvadj"
data.frame(pmvadj=pmvadj)
}
