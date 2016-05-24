# Functions return:
#
# et* - New Effective Temperature
# set* - standard Effective Temperature
# tsens - Predicted Thermal sensation 
# disc - Predicted discomfort
# ps - Predicted Percent satisfied with the Level of Air movement
# pd - Predicted Percent Dissatisfied Due to Draft
#
# File contains 1 function:
#   - calc2Node (ta, vel, rh, clo=.5, met=1, wme=0, pb=760, ltime=60, ht=166, wt=70)
#       returns et, set, tsens, disc, ps, pd
#
# The code for calc2Node is based on the code in BASIC and C++ presented by Fountain and Huizenga. 

# Function: Two Node model ################
###########################################
calc2Node <- function(ta, tr, vel, rh, clo = .5, met = 1, wme = 0, pb = 760, ltime = 60, ht = 171, wt = 70, tu = 40, obj = "set", csw = 170, cdil = 120, cstr = .5){
   
  m <- met * 58.2 #[w/m2]
  w <- wme * 58.2 #[w/m2]
  mw <- m - w 
  
  kclo <- .25	
  csw  <- 170   #driving coefficient for regulatory sweating
  cdil <- 120   #driving coefficient for vasodilation
  cstr <- .5    #driving coefficient for vasoconstriction
  
  tskn  <- 33.7  #setpoint (neutral) value for tsk
  tcrn  <- 36.8  #setpoint value for tcr
  tbn   <- 36.49 #setpoint for tb (.1*tskn + .9*tcrn)
  skbfn <- 6.3   #neutral value for skbf
  sbc   <- 5.6697 * 10 ^ (-08) #stephan-Boltzmann constant
  sa    <- ((ht * wt) / 3600 ) ^ .5   # surface Area (m2) according to mosteller formula 
  
  if (vel < .1){vel <- .1}
  
  # INITIAL VALUEs - start of 1st experiment
  tsk   <- tskn
  tcr   <- tcrn
  skbf  <- skbfn
  mshiv <- 0
  alfa  <- .1
  rmm    <- m
  esk   <- .1 * met
  
  # UNIT CONVERSIONS (from input variables)
  atm   <- pb / 760 # input unit is torr!
  timeh <- ltime / 60 
  rcl   <- .155 * clo 
  facl  <- 1 + .15 * clo  # INCreasE IN BODY sURFAce Area DUE to cloTHING
  lr    <- 2.2 / atm     # Lewis Relation is 2.2 at sea level
  pa    <- rh * exp(18.6686 - (4030.183/(ta + 235))) / 100
  
  if (clo <= 0){
    wcrit <- .38 * vel ^ (-.29)
    icl   <- 1 
  } else {
    wcrit <- .59 * vel ^ (-.08)
    icl   <- .45
  }
  
  chc <- 3 * atm  ^ .53
  if(obj == "set"){ # not used for calculation of adjusted pmv
	  if(rmm/58.2 <.85){
		chca<-0
	  } else {
		chca<-5.66*((rmm/58.2-.85)*atm) ^ .39
	  }
  }
  chcv <- 8.600001 * (vel * atm) ^ .53
  if(obj == "set"){ # not used for calculation of adjusted pmv
	if(chc<=chca){chc<-chca}
  }
  chc <- max(chc, chcv) 
  
  # initial estimate of tcl
  chr <- 4.7
  ctc <- chr + chc
  ra <- 1 / (facl * ctc) #resistance of air layer to dry heat transfer
  top <- (chr * tr + chc * ta) / ctc
  tcl <- top + (tsk - top) / (ctc * (ra + rcl))
  tim <- 1
  
  ################# BEGIN ITERATION
  
  # tcl and chr are solved iteritively using: H(tsk - to) <- ctc(tcl - to), 
  # where H <- 1/(ra + rcl) and ra <- 1/facl*ctc
  
  fnsvp  <- function(T){exp(18.6686 - 4030.183 / (T + 235))}
  #fnp    <- function(x){x * (-1) * (x > 0)}
  fnfar  <- function(T){9 / 5 * T + 32}
  fnerre <- function(x, hsk, hd, tsk, w, he, pssk){hsk - hd * (tsk - x) - w * he * (pssk - .5 * fnsvp(x))}
  fnerrs <- function(x, hsk, hd.s, tsk, w, he.s, pssk){hsk - hd.s * (tsk - x) - w * he.s * (pssk - .5 * fnsvp(x))}
  
  tclold <- tcl
  flag   <- FALSE
  for (tim in 1:ltime){
    if (flag){
      tcl <- (ra * tsk + rcl * top) / (ra + rcl)
      if (abs(tcl - tclold) > .01){
        flag   <- FALSE
        tclold <- tcl
      } else {
        flag <- TRUE
      }
    }
    while (!flag){
      chr <- 4 * sbc * ((tcl + tr) / 2 + 273.15) ^ 3 * .72
      ctc <- chr + chc
      ra  <- 1 / (facl * ctc) # resistance of air layer to dry heat transfer
      top <- (chr * tr + chc * ta) / ctc
      tcl <- (ra * tsk + rcl * top) / (ra + rcl)
      if (abs(tcl - tclold) > .01){
        flag   <- FALSE
        tclold <- tcl
      } else {
        flag   <- TRUE
      }
    } #IF ABs(tcl-tclold)>.01 TheN tcl.OLD<-tcl: GOto 1160
    dry  <- (tsk - top) / (ra + rcl)
    hfcs <- (tcr - tsk) * (5.28 + 1.163 * skbf)
    eres <- .0023 * m * (44 - pa)
    cres <- .0014 * m * (34 - ta)
    scr  <- m - hfcs - eres - cres - w # m and w in w/2
    ssk  <- hfcs - dry - esk
    tcsk <- .97 * alfa * wt
    tccr <- .97 * (1 - alfa) * wt
    dtsk <- (ssk * sa) / tcsk / 60 # deg C per minute
    dtcr <- scr * sa / tccr / 60 # deg C per minute
    dtim <- 1 #minutes
    # U<-ABs(dtsk): IF U>.1 TheN dtim<-.1/U
    # U<-ABs(dtcr): IF U>.1 AND .1/U<dtim TheN dtim<-.1/U
    tsk <- tsk + dtsk * dtim
    tcr <- tcr + dtcr * dtim
    tb  <- alfa * tsk + (1 - alfa) * tcr
    #sksig <- tsk - tskn
    #warms<-fnp(sksig)
    #colds<-fnp(-sksig)
    if (tsk > tskn){
      warms <- tsk - tskn
      colds <- 0
    } else {
      colds <- tskn - tsk
      warms <- 0
    }
    #crsig <- (tcr - tcrn)
    #warmc<-fnp(crsig)
    #coldc<-fnp(-crsig)
    if (tcr > tcrn){
      warmc <- tcr - tcrn
      coldc <- 0
    } else {
      coldc <- tcrn - tcr
      warmc <- 0
    }
    #bdsig<-tb-tbn
    #warmb<-fnp(bdsig)
    #coldb<-fnp(-bdsig)
    if (tb > tbn){
      warmb <- tb - tbn
      coldb <- 0
    } else {
      coldb <- tbn - tb
      warmb <- 0
    }
    skbf <- (skbfn + cdil * warmc) / (1 + cstr * colds)
    if (skbf > 90){skbf <- 90}
    if (skbf < .5){skbf <- .5}
    regsw <- csw * warmb * exp(warms / 10.7)
    if (regsw > 500){regsw <- 500}
    ersw <- .68 * regsw
    # lr <- 2.02*(tsk+273.15)/273.15
    rea  <- 1 / (lr * facl * chc) #evaporative resistance of air layer
    recl <- rcl / (lr * icl) #evaporative resistance of clothing (icl<-.45)
    emax <- (fnsvp(tsk) - pa) / (rea + recl)
    prsw <- ersw / emax
    pwet <- .06 + .94 * prsw
    edif <- pwet * emax - ersw
    esk  <-ersw + edif
    if (pwet > wcrit){
      pwet <- wcrit
      prsw <- (wcrit) / .94 # (wcrit-.06)/.94
      ersw <- prsw * emax
      edif <- .06 * (1 - prsw) * emax
      esk  <- ersw + edif
    }
    if (emax < 0){
      edif <- 0
      ersw <- 0
      pwet <- wcrit
      prsw <- wcrit
      esk  <- emax
    }
    esk   <- ersw + edif
    mshiv <- 19.4 * colds * coldc
    m     <- rmm + mshiv
    alfa  <- .0417737 + .7451833 / (skbf + .585417)
    #GOsUB2680 'screen output
    #IF OUtopT<-1 TheN GOsUB 3900 'minute-by-minute hcopy output
    #tim <- tim+dtim #IF tim<=ltime TheN GOto 1200
  }
  
  #####################################################################
  # CALCULATE COMFORT INDICES
  #####################################################################
  
  # Define new heat flow terms, coeffs, and abbreviations
  store   <- m - w - cres - eres - dry - esk     #rate of body heat storage #?
  hsk     <- dry + esk                   #total heat loss from skin
  rn      <- m - w                       #net metabolic heat production [w/m2]
  ecomf   <- .42 * (rn - 58.2)
  
  if (ecomf < 0){ecomf <- 0} #from Fanger
  
  ereq    <- rn - eres - cres - dry#?
  emax    <- emax * wcrit
  hd      <- 1 / (ra + rcl) #?
  he      <- 1 / (rea + recl)#?
  wet     <- pwet
  pssk    <- fnsvp(tsk)
  
  #Definition of ASHRAE standard environment... denoted "s"
  chrs <- chr
  
  if (met < .85){
    chcs <- 3
  } else { # GOto 1950
    chcs <- 5.66 * ((met - .85)) ^ .39
    if (chcs < 3){chcs <- 3}
  }
  ctcs  <- chcs + chrs
  rclos <- 1.52 / ((met - wme) + .6944) -.1835 # here wme als in [met]!
  rcls  <- .155 * rclos
  facls <- 1 + kclo * rclos
  fcls  <- 1 / (1 + .155 * facls * ctcs * rclos)
  ims   <- .45
  icls  <- ims * chcs / ctcs * (1 - fcls) / (chcs / ctcs - fcls * ims)
  ras   <- 1 / (facls * ctcs)
  reas  <- 1 / (lr * facls * chcs)
  recls <- rcls / (lr * icls)
  hd.s  <- 1 / (ras + rcls)
  he.s  <- 1 / (reas + recls)
  
  # et* (standardized humidity/ actual do, pb, and chc)
  # determined using Newton's iterative solution
  # fnerre is defined in general setup section above
  delta <- .0001
  xold <- tsk - hsk / hd #lower bound for et*
  
  flag1 <- FALSE
  while (!flag1){
    err1 <- fnerre(xold, hsk, hd, tsk, wet, he, pssk)
    err2 <- fnerre(xold + delta, hsk, hd, tsk, wet, he, pssk)
    x <- xold - delta * err1 / (err2 - err1)
    if (abs(x - xold) > .01){
      xold <- x
      flag1 <- FALSE
    } else {
      flag1 <- TRUE
    }
  } #IF ABs (x-xold)>.01 TheN xold<-x: GOto 2120
  et <- x
  
  # set* (standardized humidity, clo, pb, and chc)
  # determined using Newton's iterative solution
  # fnerrs is defined in the GENEraL setuP section above
  
  xold <- tsk-hsk / hd.s #lower bound for set
  flag2 <- FALSE
  
  while (!flag2){
    err1 <- fnerrs(xold, hsk, hd.s, tsk, wet, he.s, pssk)
    err2 <- fnerrs(xold + delta, hsk, hd.s, tsk, wet, he.s, pssk)
    
	x <- xold - delta * err1 / (err2 - err1)
    
	if (abs(x - xold) > .01){
      xold <- x
      flag2 <- FALSE
    } else {
      flag2 <- TRUE
    }
  } # IF ABs(x-xold)>.01 TheN xold<-x: GOto 2220
  set <- x
  
  # sto <- standard operative temperature:
  # defined by: (tsk-sto)/(ras+rcls)<-(tsk-to)/(ra + rcl)
  sto <- tsk - (ras + rcls) * (tsk - top) / (ra + rcl)
  
  # sVPO <- standard operative vapor pressure:
  # defined by: (psk - sVPO)/(reas+recls) <- (psk-pa)/(rea+recl)
  # sVPO <- psK - (reas+recls)*(psK-pa)/(rea+recl)
  
  # tsens is a function of tb 
  tbml <- (.194 / 58.15) * rn + 36.301 #lower limit for evaporative regulation
  tbmh <- (.347 / 58.15) * rn + 36.669 #upper limit for evaporative regulation
  if (tb < tbml){
    tsens <- .4685 * (tb - tbml)
  } else if (tb >= tbml & tb < tbmh){
    tsens <- wcrit * 4.7 * (tb - tbml) / (tbmh - tbml)
  } else if (tb >= tbmh){
    tsens <- wcrit * 4.7 + .4685 * (tb - tbmh)
  }
  
  # disc varies with relative thermoregulatory heat strain
  # Valid only when disc>0. when disc<0, disc is numerically <-tsens.
  disc <- 4.7 * (ersw - ecomf) / (emax - ecomf - edif)
  if (disc < 0){disc <- tsens}
  
  # Calculate Gagge's version of Fanger's Predicted mean Vote (pmv)
  pmvg <- (.303 * exp(-.036 * m) + .028) * (ereq - ecomf - edif)
  
  #Gagge's pmv.set is the same as Fanger's pmv except that dry is calculated
  #using set* rather than top
  dry2  <- hd.s * (tsk - set)
  ereq2 <- rn - cres - eres - dry2
  pmv2  <- (.303 * exp(-.036 * m) + .028) * (ereq2 - ecomf - edif)
  
  # other indices related to air movement
  tue  <- (34 - ta) * (vel - 0.05) ^ 0.6223
  pd   <- tue * (3.143 + 0.3696 * vel * tu)
  ps   <- 100 * (1.13 * (top ^ .5) - .24 * top + 2.7 * (vel ^ .5) - .99 * vel)
  pts  <- .25 * set - 6.03
  
  output <- data.frame(et, set, tsens, disc, pd, ps, pts)
  rm(et, set, tsens, disc, pd, ps, pts)
  output
}

 
## calc aset ################################
## based on Gao et al.
############################################
calcPtsa <- function(ta, tr, vel, rh, clo = .5, met = 1, wme = 0, pb = 760, ltime = 60, ht = 171, wt = 70, tu = 40, asCoeff){
 
	set <- calc2Node(ta, tr, vel, rh, clo, met, wme, pb, ltime, ht, wt, tu, obj = "set")[2]
	ptsa <- .25 * set - 6.03
	ptsa <- ptsa / (1 + asCoeff * ptsa)
	names(ptsa) <- "ptsa"
	data.frame(ptsa = ptsa)

}

## calc eset ################################
## based on Gao et al.
############################################
calcPtse <- function(ta, tr, vel, rh, clo = .5, met = 1, wme = 0, pb = 760, ltime = 60, ht = 171, wt = 70, tu = 40, esCoeff, asv){

	met  <- ifelse (asv > 0, met * (1 + asv * (-.067)), met)
	set  <- calc2Node(ta, tr, vel, rh, clo, met, wme, pb, ltime, ht, wt, tu, obj = "set")[2]
	ptse <- .25 * set - 6.03
	ptse <- esCoeff * ptse
	names(ptse) <- "ptse"
	data.frame(ptse = ptse)

}


 