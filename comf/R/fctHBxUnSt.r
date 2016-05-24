## R program to calculate unsteady state human body exergy consumPtion

##    This is a program for the calculation of human-body core and skin temperatures
##    and also clothing surface temperature based on the two-node model
##    originally developed by Gagge et al.
##    The program has been developed so that it fits the calculation of human-body
##    exergy balance under unsteady-state conditions.
##    The program is based on the Excel version for calculating human body exergy consumPtion rate developed by masanori Shukuya
##                                        1st ver.  masanori Shukuya 13th February, 2013
##
##    This program has been further extended to be able to include the human-body exergy balance.
##                                                  masanori Shukuya 11th may, 2014
##
##    This version is for un-steady state exergy calculation.
##                                                  masanori Shukuya 30th June, 2014//18th February, 2015
##
## 	 transformation of VBA-code and Excel procedures into R syntax
## 													Marcel Schweiker May, 2015

## comment to values of cdil and sigmatr: According to Gagges paper (1973), the value of cdil may vary between 75 and 225 and sigma-tr between 0.25 and 0.75. There is a note in the appendix of his paper saying two things:  1) whatever the values taken for cdil and sigma-tr, there must be no significant change in resulting thermal equilibrium. But, the values taken for cdil and sigmatr do affect time to equilibrium. According to our analysis, the values of 100 and .25 lead to the best fit of calculated and observed skin temperature. Therefore, this value was used for the current paper.

#############################################
## Main program
#############################################

## This is a program for the calculation of human-body core and skin temperatures
## and also clothing surface temperature based on the two-node model
## originally developed by Gagge et al.
## The program has been developed so that it fits the calculation of human-body
## exergy balance under unsteady-state conditions.
## 1st ver.  ms 13th February, 2013

## This program has been further extended to be able to include the human-body exergy balance.
##  ms 11th may, 2014

## This version is for un-steady state exergy calculation.
##  ms 30th June, 2014//18th February, 2015

## transfer of VBA-code to R-code by m. schweiker march, 2015

###############################################

calcHbExUnsteady <- function(ta, tr, rh, vel, clo, met, tao, rho, frad = .7, eps = .95,  ic = 1.085, ht=171, wt=70, tcr=37, tsk=36, basMet=58.2, warmUp=60, cdil=100, sigmatr=.25, dateTime){

## definition of output variables
		## Exergy input
		xInmetu <- xInmetwcu <- xInAIRwcu <- xInAIRwcwcu <- xInAIRwdu <- xInAIRwdwdu <- xInLUNGwcu <- xInLUNGwcwcu <- xInLUNGwdu <- xInLUNGwdwdu <- xInsheLLwcu <- xInsheLLwcwcu <- xInsheLLwdu <- xInsheLLwdwdu <- xInraDu <- xInraDwcu <- xIntotaLu <- NA
		
		## Exergy output
		xoutstorecoreu <- xoutstoreshelu <- xoutaIRwcu <- xoutaIRwcwcu <- xoutaIRwdu <- xoutaIRwdwdu <- xoutswEATwcu <- xoutswEATwcwcu <- xoutswEATwdu <- xoutswEATwdwdu <- xoutraDu <- xoutraDwcu <- xoutCONVu <- xoutCONVwcu <- xouttotaLu <- NA
		
		## balance and additional variables
		xconsu <- tsku <- tcru <- wu <- NA

	tko <- 273.15
	tskSet <- 33.7
	tcrSet <- 36.8
	lr <- 16.5 * 10  ^  (-3)
	csw <- 170
	mAir <- 28.97 * 0.001; 
	mWater <- 18.015 * 0.001; 
	rGas <- 8.31446 #[J/(molK)]
	Row <- 1000; Roa <- 1.2#[kg/m3]
	cpBody <- 3490; cpa <- 1005; cpv <- 1846; cpw <- 4186 #[J/(kgK)]
	PO <- 101325#[N/m2<-J/m3]
	
	## radiation area of human body - taking account for covered parts by other parts of body (e.g. inner parts of arm) - coefficients taken from Fanger
	aBody <- 0.008883 * ht  ^  0.663 * wt  ^  0.444
	rMA <- wt / aBody 
		
	## A series of calculation timestep by timestep
	##  adjust here for input of individual conditions in each line

	runs <- length(ta)

	i <- 1
		
		icl <- clo[i] # clothing insulation [clo]
		va <- vel[i] # Indoor air velocity [m/s]
		metrate <- met[i] # metabolic rate [met]
		tmr <- tr[i] # mean radiant temp [degree C] 
		tair <- ta[i] # air temp [degree C] 
		phia <- rh[i] # relative humidity indoors [%]
		toEnv <- tao[i] # outdoor temp [degree C] 
		phiaEnv <- rho[i] # relative humidity outdoors [%]

	## warm-up period - see also discussion in paper	
	for (j in 1:warmUp){
		# extract single value from time series data
		dtCal <- 60 
		qmet <- metaTherm(metrate, basMet)
		fcl <- 1 + 0.15 * icl
		rcl <- 0.155 * icl

		TKa <- tair + tko; tkmr <- tmr + tko

		tkoEnv <- toEnv + tko 
		pvEnv <- pVapor(toEnv, phiaEnv) # function see Book Shukuya p. 268

		hr <- 6.13 * frad * eps # radiative heat transfer coefficient
		hc <- hcvG(va, metrate, basMet) # convective heat transfer coefficient
		top <- (hr * tmr + hc * tair) / (hr + hc) 
		ffcl <- 1 / (1 + 0.155 * icl * fcl * (hr + hc)); 
		fpcl <- 1 / (1 + 0.155 * icl * fcl * hc / ic) # related to evaporation
		cl <- 1 / (1 + 0.155 * icl * fcl * hc / ic)
		im <- hc * fpcl / ((hr + hc) * ffcl)

		iclStar <- 0.6; fclStar <- 1 + 0.15 * iclStar
		hcStar <- hcvG(0.1, 1, basMet)
		ffclStar <- 1 / (1 + 0.155 * iclStar * fclStar * (hr + hcStar))
		fpclStar <- 1 / (1 + 0.155 * iclStar * fclStar * hcStar / ic)
		imStar <- hcStar * fpclStar / ((hr + hcStar) * ffclStar); 
		cres <- 0.0014 * (basMet * metrate) * (34 - tair) # heat transfer to environment by convection in relation to respiration (empirical equation by Fanger)

		pva <- pVapor(tair, phia)
		eres <- 0.0173 * (basMet * metrate) * (5.87 - pva / 1000) # heat removed by evaporation in relation to respiration (empirical equation by Fanger)

		
		qShiv <- mshiv(tcrSet, tskSet, tcr, tsk)
		vblS <- vblCdilStr(cdil, sigmatr, tcrSet, tskSet, tcr, tsk)# blood flow rate
		# next line is core of Gagge model. If blood flow becomes lower, than skin layer gets more dominant
		alfaSk <- 0.0418 + 0.745 / (vblS + 0.585); # factor to adjust for difference in dominance
		kS <- 5.28 + 1.163 * vblS # conductance between core and skin layer
		Qcr <- (1 - alfaSk) * rMA * cpBody; Qsk <- alfaSk * rMA * cpBody # heat capacity of core and skin layer
		
		tcrN <- (1 - dtCal * kS / Qcr) * tcr + dtCal / Qcr * (qmet + qShiv - (cres + eres) + kS * tsk) # core temperature at time step before step calculated
		psks <- pVapor(tsk, 100)# saturated water vapour pressure at skin surface calculated with pure water , i.e. adaptive processes such as less salt in sweat might be put heresee also p. 281 of book Shukuya
		emax <- fpcl * (fcl * lr * hc) * (psks - pva) # max rate of water dispersion when the whole skin surface is wet
		tb <- alfaSk * tsk + (1 - alfaSk) * tcr; # average body temperature using calculated tsk and tcr
		tbSet <- alfaSk * tsk + (1 - alfaSk) * tcrSet # average body temperature using set-point values for tsk and tcr
		ersw <- csw * (tb - tbSet) * exp((tsk - tskSet) / 10.7) # amount of sweat secretion 
		
		w <- 0.06 + 0.94 * ersw / emax
		if (w < 0.06){
			w <- 0.06
		}
		if (1 < w){
			w <- 1
		}
		
		DTQ <- dtCal / Qsk
		tskN <- (1 - DTQ * kS - DTQ * fcl * ffcl * (hr + hc)) * tsk - DTQ * w * fcl * lr * hc * fpcl * psks + DTQ * (kS * tcr + fcl * (hr + hc) * ffcl * top + w * fcl * lr * hc * fpcl * pva) # tsk in next step
		tclN <- ((1 / rcl) * tskN + fcl * hr * tmr + fcl * hc * tair) / (1 / rcl + fcl * hr + fcl * hc)
		
		## Exergy balance
		
		## Thermal exergy generation by metabolism
		
		tkcrN <- tcrN + tko; tkskN <- tskN + tko; tkclN <- tclN + tko
		metaEnergy <- qmet + qShiv
		xMet <- metaEnergy * (1 - tkoEnv / tkcrN); 
		xMetwc <- wcXCheck(tkcrN, tkoEnv)

		## Inhaled humid air

		Vin <- 1.2 * 10  ^  (-6) * metaEnergy # Volume of air intake [V/s]
		cpav <- cpa * (mAir / (rGas * TKa)) * (PO - pva) + cpv * (mWater / (rGas * TKa)) * pva
		xInhaleWc <- Vin * wcEx(cpav, TKa, tkoEnv); 
		xInhaleWcwc <- wcXCheck(TKa, tkoEnv)

		xInhaleWd <- Vin * wdEx(TKa, tkoEnv, pva, pvEnv); 
		xInhaleWdwd <- wdXCheck(pva, pvEnv)
		
		## Liquid water generated in the core by metabolism to be dispersed into the exhaled air
		
		VwCore <- Vin * (0.029 - 0.049 * 10  ^  (-4) * pva)
		xLwCoreWc <- VwCore * Roa * wcEx(cpw, tkcrN, tkoEnv); xLwCoreWcwc <- wcXCheck(tkcrN, tkoEnv)
		
		pvs_env <- pVapor(toEnv, 100)
		xLwCoreWet <- VwCore * Roa * rGas / mWater * tkoEnv * log(pvs_env / pvEnv)
		xLwCoreWetwd <- "wet"
		
		## Liquid water generated in the shell by metabolism to be secreted as sweat
		
		vwShellRow <- w * emax / (2450 * 1000)
		xLwShellWc <- vwShellRow * wcEx(cpw, tkskN, tkoEnv); xLwShellWcwc <- wcXCheck(tkskN, tkoEnv)
		
		xLwShellWd <- vwShellRow * wdExLw(tkoEnv, pvs_env, pva, pvEnv); 
		xLwShellWdwd <- wdXCheck(pva, pvEnv)
		
		## radiant exergy input

		xInRad <- fcl * hr * (tkmr - tkoEnv)  ^  2 / (tkmr + tkoEnv); 
		xInRadwc <- wcXCheck(tkmr, tkoEnv)

		## total exergy input

		xIntotal <- xMet + xInhaleWc + xInhaleWd + xLwCoreWc + xLwCoreWet + xLwShellWc + xLwShellWd + xInRad

		## Exergy stored

		xStCore <- Qcr * (1 - tkoEnv / tkcrN) * (tcrN - tcr) / dtCal
		xStShell <- Qsk * (1 - tkoEnv / tkskN) * (tskN - tsk) / dtCal

		## Exhaled humid air

		pvssCr <- pVapor(tcrN, 100)
		cpav <- cpa * (mAir / (rGas * tkcrN)) * (PO - pvssCr) + cpv * (mWater / (rGas * tkcrN)) * pvssCr
		xExhaleWc <- Vin * wcEx(cpav, tkcrN, tkoEnv); 
		xExhaleWcwc <- wcXCheck(tkcrN, tkoEnv) 

		xExhaleWd <- Vin * wdEx(tkcrN, tkoEnv, pvssCr, pvEnv); xExhaleWdwd <- wdXCheck(pvssCr, pvEnv)

		## water vapor originating from the sweat and humid air containing the evaporated sweat

		xSweatWc <- vwShellRow * wcEx(cpv, tkclN, tkoEnv); 
		xSweatWcwc <- wcXCheck(tkclN, tkoEnv)

		xSweatWd <- vwShellRow * wdExLw(tkoEnv, pva, pva, pvEnv); 
		xSweatWdwd <- wdXCheck(pva, pvEnv) 

		## radiant exergy output

		xOutRad <- fcl * hr * (tkclN - tkoEnv)  ^  2 / (tkclN + tkoEnv); 
		xOutRadwc <- wcXCheck(tkclN, tkoEnv)

		## Exergy transfer by convection

		xOutConv <- fcl * hc * (tkclN - TKa) * (1. - tkoEnv / tkclN); 
		xOutConvwc <- wcXCheck(tkclN, tkoEnv)

		xouttotal <- xStCore + xStShell + xExhaleWc + xExhaleWd + xSweatWc + xSweatWd + xOutRad + xOutConv

		xConsumption <- xIntotal - xouttotal

		#etStar <- calcet(top, tair, phia, w, im, 50, im);
		
		tcr <- tcrN
		tsk <- tskN	
	
	}
	
		## Output values
		## Exergy input
		xInmetu[i] <- xMet #metabolism 
		xInmetwcu[i] <- xMetwc #metabolism warm/cold
		xInAIRwcu[i] <- xInhaleWc	# inhaled humid air
		xInAIRwcwcu[i] <- xInhaleWcwc
		xInAIRwdu[i] <- xInhaleWd # 
		xInAIRwdwdu[i] <- xInhaleWdwd # wet/dry
		xInLUNGwcu[i] <- xLwCoreWc # water lung
		xInLUNGwcwcu[i] <- xLwCoreWcwc
		xInLUNGwdu[i] <- xLwCoreWet
		xInLUNGwdwdu[i] <- xLwCoreWetwd
		xInsheLLwcu[i] <- xLwShellWc # water sweat
		xInsheLLwcwcu[i] <- xLwShellWcwc
		xInsheLLwdu[i] <- xLwShellWd
		xInsheLLwdwdu[i] <- xLwShellWdwd
		xInraDu[i] <- xInRad # radiation in
		xInraDwcu[i] <- xInRadwc
		xIntotaLu[i] <- xIntotal # totaL exergy input
		
		## Exergy output
		xoutstorecoreu[i] <- xStCore # exergy stored
		xoutstoreshelu[i] <- xStShell
		xoutaIRwcu[i] <- xExhaleWc # exhaled humid air
		xoutaIRwcwcu[i] <- xExhaleWcwc
		xoutaIRwdu[i] <- xExhaleWd
		xoutaIRwdwdu[i] <- xExhaleWdwd
		xoutswEATwcu[i] <- xSweatWc # water vapour
		xoutswEATwcwcu[i] <- xSweatWcwc
		xoutswEATwdu[i] <- xSweatWd 
		xoutswEATwdwdu[i] <- xSweatWdwd
		xoutraDu[i] <- xOutRad # radiation out
		xoutraDwcu[i] <- xOutRadwc
		xoutCONVu[i] <- xOutConv # convection
		xoutCONVwcu[i] <- xOutConvwc
		xouttotaLu[i] <- xouttotal # totaL exergy out
		
		## balance
		xconsu[i] <- xConsumption # exergy consumPtion total
		tsku[i] <- tsk
		tcru[i] <- tcr
		wu[i] <- w
	
	## from here start dynamic series of data from row 2 of dataset
	for (i in 2:runs){ 
	
		## Calc Dtcal dependend on timediff between rows + approx if more than 20 seconds time difference
		dtCal <-  as.numeric(difftime(dateTime[i], dateTime[i-1], units="secs")) # in seconds 
		
		if(dtCal > 20){ # if time difference is greater than 20 seconds (leading to the calcs getting unstable), intermediate values or calculated
		
			nRunsInter <- trunc(dtCal/20-1, 0) # calc number of minutes in between interval
			dtCal <- dtCal / (nRunsInter+1)
			
			for (k in 1: nRunsInter) {
				
				## get linear approximated values for intermediate times
				icl <- clo[i-1] + k*((clo[i] - clo[i-1])/(nRunsInter+1)) # clothing insulation [clo]
				va <- vel[i-1] + k*((vel[i] - vel[i-1])/(nRunsInter+1)) # Indoor air velocity [m/s]
				metrate <- met[i-1] + k*((met[i] - met[i-1])/(nRunsInter+1)) # metabolic rate [met]
				tmr <- tr[i-1] + k*((tr[i] - tr[i-1])/(nRunsInter+1)) # mean radiant temp [degree C] 
				tair <- ta[i-1] + k*((ta[i] - ta[i-1])/(nRunsInter+1)) # air temp [degree C] 
				phia <- rh[i-1] + k*((rh[i] - rh[i-1])/(nRunsInter+1)) # relative humidity indoors [%]
				toEnv <- tao[i-1] + k*((tao[i] - tao[i-1])/(nRunsInter+1)) # outdoor temp [degree C] 
				phiaEnv <- rho[i-1] + k*((rho[i] - rho[i-1])/(nRunsInter+1)) # relative humidity outdoors [%]
				
				qmet <- metaTherm(metrate, basMet)
				fcl <- 1 + 0.15 * icl
				rcl <- 0.155 * icl

				TKa <- tair + tko; tkmr <- tmr + tko

				tkoEnv <- toEnv + tko 
				pvEnv <- pVapor(toEnv, phiaEnv) #--------function see Book Shukuya p. 268

				hr <- 6.13 * frad * eps # radiative heat transfer coefficient
				hc <- hcvG(va, metrate, basMet) # convective heat transfer coefficient 
				top <- (hr * tmr + hc * tair) / (hr + hc) 
				ffcl <- 1 / (1 + 0.155 * icl * fcl * (hr + hc)); 
				fpcl <- 1 / (1 + 0.155 * icl * fcl * hc / ic) 
				cl <- 1 / (1 + 0.155 * icl * fcl * hc / ic)
				im <- hc * fpcl / ((hr + hc) * ffcl)
				
				iclStar <- 0.6; fclStar <- 1 + 0.15 * iclStar
				hcStar <- hcvG(0.1, 1, basMet)
				ffclStar <- 1 / (1 + 0.155 * iclStar * fclStar * (hr + hcStar))
				fpclStar <- 1 / (1 + 0.155 * iclStar * fclStar * hcStar / ic)
				imStar <- hcStar * fpclStar / ((hr + hcStar) * ffclStar); 
				cres <- 0.0014 * (basMet * metrate) * (34 - tair) # heat transfer to environment by convection in relation to respiration (empirical equation by Fanger)

				pva <- pVapor(tair, phia)
				eres <- 0.0173 * (basMet * metrate) * (5.87 - pva / 1000) # heat removed by evaporation in relation to respiration (empirical equation by Fanger)

				qShiv <- mshiv(tcrSet, tskSet, tcr, tsk)
				vblS <- vblCdilStr(cdil, sigmatr, tcrSet, tskSet, tcr, tsk)# blood flow rate 
				## next line is core of Gagge model. If blood flow becomes lower, than skin layer gets more dominant
				alfaSk <- 0.0418 + 0.745 / (vblS + 0.585); # factor to adjust for difference in dominance
				kS <- 5.28 + 1.163 * vblS # conductance between core and skin layer
				Qcr <- (1 - alfaSk) * rMA * cpBody; Qsk <- alfaSk * rMA * cpBody # heat capacity of core and skin layer
				
				tcrN <- (1 - dtCal * kS / Qcr) * tcr + dtCal / Qcr * (qmet + qShiv - (cres + eres) + kS * tsk) # core temperature at time step before step calculated
				psks <- pVapor(tsk, 100)# saturated water vapour pressure at skin surface calculated with pure water , i.e. adaptive processes such as less salt in sweat might be put heresee also p. 281 of book Shukuya
				emax <- fpcl * (fcl * lr * hc) * (psks - pva) # max rate of water dispersion when the whole skin surface is wet
				tb <- alfaSk * tsk + (1 - alfaSk) * tcr; # average body temperature using calculated tsk and tcr
				tbSet <- alfaSk * tsk + (1 - alfaSk) * tcrSet # average body temperature using set-point values for tsk and tcr
				ersw <- csw * (tb - tbSet) * exp((tsk - tskSet) / 10.7) # amount of sweat secretion 
				
				w <- 0.06 + 0.94 * ersw / emax
				if (w < 0.06){
					w <- 0.06
				}
				if (1 < w){
					w <- 1
				}
				
				DTQ <- dtCal / Qsk
				tskN <- (1 - DTQ * kS - DTQ * fcl * ffcl * (hr + hc)) * tsk - DTQ * w * fcl * lr * hc * fpcl * psks + DTQ * (kS * tcr + fcl * (hr + hc) * ffcl * top + w * fcl * lr * hc * fpcl * pva) # tsk in next step
				tclN <- ((1 / rcl) * tskN + fcl * hr * tmr + fcl * hc * tair) / (1 / rcl + fcl * hr + fcl * hc)
				
				#etStar <- calcet(top, tair, phia, w, im, 50, im); 
				
				tcr <- tcrN
				tsk <- tskN	
					
			} # end for k times in between interval
		
		} # end if Dtcal >20

		## extract single value from time series data
		icl <- clo[i] # clothing insulation [clo]
		va <- vel[i] # Indoor air velocity [m/s]
		metrate <- met[i] # metabolic rate [met]
		tmr <- tr[i] # mean radiant temp [degree C] 
		tair <- ta[i] # air temp [degree C] 
		phia <- rh[i] # relative humidity indoors [%]
		toEnv <- tao[i] # outdoor temp [degree C] 
		phiaEnv <- rho[i] # relative humidity outdoors [%]
			
		qmet <- metaTherm(metrate, basMet)
		fcl <- 1 + 0.15 * icl
		rcl <- 0.155 * icl

		TKa <- tair + tko; tkmr <- tmr + tko

		tkoEnv <- toEnv + tko 
		pvEnv <- pVapor(toEnv, phiaEnv) # function see Book Shukuya p. 268

		hr <- 6.13 * frad * eps # radiative heat transfer coefficient
		hc <- hcvG(va, metrate, basMet) # convective heat transfer coefficient #
		top <- (hr * tmr + hc * tair) / (hr + hc) 
		ffcl <- 1 / (1 + 0.155 * icl * fcl * (hr + hc)); 
		fpcl <- 1 / (1 + 0.155 * icl * fcl * hc / ic) # related to evaporation
		cl <- 1 / (1 + 0.155 * icl * fcl * hc / ic)
		im <- hc * fpcl / ((hr + hc) * ffcl)
		
		iclStar <- 0.6; fclStar <- 1 + 0.15 * iclStar
		hcStar <- hcvG(0.1, 1, basMet)
		ffclStar <- 1 / (1 + 0.155 * iclStar * fclStar * (hr + hcStar))
		fpclStar <- 1 / (1 + 0.155 * iclStar * fclStar * hcStar / ic)
		imStar <- hcStar * fpclStar / ((hr + hcStar) * ffclStar); 
		cres <- 0.0014 * (basMet * metrate) * (34 - tair) # heat transfer to environment by convection in relation to respiration (empirical equation by Fanger)

		pva <- pVapor(tair, phia)
		eres <- 0.0173 * (basMet * metrate) * (5.87 - pva / 1000) # heat removed by evaporation in relation to respiration (empirical equation by Fanger)

		qShiv <- mshiv(tcrSet, tskSet, tcr, tsk)
		vblS <- vblCdilStr(cdil, sigmatr, tcrSet, tskSet, tcr, tsk)# blood flow rate
		# next line is core of Gagge model. If blood flow becomes lower, than skin layer gets more dominant
		alfaSk <- 0.0418 + 0.745 / (vblS + 0.585); # factor to adjust for difference in dominance
		kS <- 5.28 + 1.163 * vblS # conductance between core and skin layer
		Qcr <- (1 - alfaSk) * rMA * cpBody; Qsk <- alfaSk * rMA * cpBody # heat capacity of core and skin layer
		
		tcrN <- (1 - dtCal * kS / Qcr) * tcr + dtCal / Qcr * (qmet + qShiv - (cres + eres) + kS * tsk) # core temperature at time step before step calculated
		psks <- pVapor(tsk, 100)# saturated water vapour pressure at skin surface calculated with pure water , i.e. adaptive processes such as less salt in sweat might be put heresee also p. 281 of book Shukuya
		emax <- fpcl * (fcl * lr * hc) * (psks - pva) # max rate of water dispersion when the whole skin surface is wet
		tb <- alfaSk * tsk + (1 - alfaSk) * tcr; # average body temperature using calculated tsk and tcr
		tbSet <- alfaSk * tsk + (1 - alfaSk) * tcrSet # average body temperature using set-point values for tsk and tcr
		ersw <- csw * (tb - tbSet) * exp((tsk - tskSet) / 10.7) # amount of sweat secretion 
				
		w <- 0.06 + 0.94 * ersw / emax
		if (w < 0.06){
			w <- 0.06
		}
		if (1 < w){
			w <- 1
		}
		
		DTQ <- dtCal / Qsk
		tskN <- (1 - DTQ * kS - DTQ * fcl * ffcl * (hr + hc)) * tsk - DTQ * w * fcl * lr * hc * fpcl * psks + DTQ * (kS * tcr + fcl * (hr + hc) * ffcl * top + w * fcl * lr * hc * fpcl * pva) # tsk in next step
		tclN <- ((1 / rcl) * tskN + fcl * hr * tmr + fcl * hc * tair) / (1 / rcl + fcl * hr + fcl * hc)
		
		## Exergy balance
		
		## Thermal exergy generation by metabolism
		
		tkcrN <- tcrN + tko; tkskN <- tskN + tko; tkclN <- tclN + tko
		metaEnergy <- qmet + qShiv
		xMet <- metaEnergy * (1 - tkoEnv / tkcrN); 
		xwc <- wcXCheck(tkcrN, tkoEnv)

		## Inhaled humid air

		Vin <- 1.2 * 10  ^  (-6) * metaEnergy # Volume of air intake [V/s]
		cpav <- cpa * (mAir / (rGas * TKa)) * (PO - pva) + cpv * (mWater / (rGas * TKa)) * pva
		xInhaleWc <- Vin * wcEx(cpav, TKa, tkoEnv); xwc <- wcXCheck(TKa, tkoEnv)

		xInhaleWd <- Vin * wdEx(TKa, tkoEnv, pva, pvEnv); xwd <- wdXCheck(pva, pvEnv)
		
		## Liquid water generated in the core by metabolism to be dispersed into the exhaled air
		
		VwCore <- Vin * (0.029 - 0.049 * 10  ^  (-4) * pva)
		xLwCoreWc <- VwCore * Roa * wcEx(cpw, tkcrN, tkoEnv); xwc <- wcXCheck(tkcrN, tkoEnv)
		
		pvs_env <- pVapor(toEnv, 100)
		xLwCoreWet <- VwCore * Roa * rGas / mWater * tkoEnv * log(pvs_env / pvEnv)
		
		## Liquid water generated in the shell by metabolism to be secreted as sweat
		
		vwShellRow <- w * emax / (2450 * 1000)
		xLwShellWc <- vwShellRow * wcEx(cpw, tkskN, tkoEnv); xwc <- wcXCheck(tkskN, tkoEnv)
		
		xLwShellWd <- vwShellRow * wdExLw(tkoEnv, pvs_env, pva, pvEnv); xwd <- wdXCheck(pva, pvEnv) 
		
		## radiant exergy input

		xInRad <- fcl * hr * (tkmr - tkoEnv)  ^  2 / (tkmr + tkoEnv); xwc <- wcXCheck(tkmr, tkoEnv)

		## total exergy input

		xIntotal <- xMet + xInhaleWc + xInhaleWd + xLwCoreWc + xLwCoreWet + xLwShellWc + xLwShellWd + xInRad

		## Exergy stored

		xStCore <- Qcr * (1 - tkoEnv / tkcrN) * (tcrN - tcr) / dtCal
		xStShell <- Qsk * (1 - tkoEnv / tkskN) * (tskN - tsk) / dtCal

		## Exhaled humid air

		pvssCr <- pVapor(tcrN, 100)
		cpav <- cpa * (mAir / (rGas * tkcrN)) * (PO - pvssCr) + cpv * (mWater / (rGas * tkcrN)) * pvssCr
		xExhaleWc <- Vin * wcEx(cpav, tkcrN, tkoEnv); xwc <- wcXCheck(tkcrN, tkoEnv)

		xExhaleWd <- Vin * wdEx(tkcrN, tkoEnv, pvssCr, pvEnv); xwd <- wdXCheck(pvssCr, pvEnv)

		## water vapor originating from the sweat and humid air containing the evaporated sweat

		xSweatWc <- vwShellRow * wcEx(cpv, tkclN, tkoEnv); xwc <- wcXCheck(tkclN, tkoEnv)

		xSweatWd <- vwShellRow * wdExLw(tkoEnv, pva, pva, pvEnv); xwd <- wdXCheck(pva, pvEnv)

		## radiant exergy output

		xOutRad <- fcl * hr * (tkclN - tkoEnv)  ^  2 / (tkclN + tkoEnv); xwc <- wcXCheck(tkclN, tkoEnv)

		## Exergy transfer by convection

		xOutConv <- fcl * hc * (tkclN - TKa) * (1. - tkoEnv / tkclN); xwc <- wcXCheck(tkclN, tkoEnv)

		xouttotal <- xStCore + xStShell + xExhaleWc + xExhaleWd + xSweatWc + xSweatWd + xOutRad + xOutConv

		xConsumption <- xIntotal - xouttotal

		#etStar <- calcet(top, tair, phia, w, im, 50, im);

		tcr <- tcrN
		tsk <- tskN
		
		## Output values
		## Exergy input
		xInmetu[i] <- signif(xMet, 4) #metabolism 
		xInmetwcu[i] <- xMetwc #metabolism warm/cold
		xInAIRwcu[i] <- signif(xInhaleWc, 4)	# inhaled humid air
		xInAIRwcwcu[i] <- xInhaleWcwc
		xInAIRwdu[i] <- signif(xInhaleWd, 4) # 
		xInAIRwdwdu[i] <- xInhaleWdwd # wet/dry
		xInLUNGwcu[i] <- signif(xLwCoreWc, 4) # water lung
		xInLUNGwcwcu[i] <- xLwCoreWcwc
		xInLUNGwdu[i] <- signif(xLwCoreWet, 4)
		xInLUNGwdwdu[i] <- xLwCoreWetwd
		xInsheLLwcu[i] <- signif(xLwShellWc, 4) # water sweat
		xInsheLLwcwcu[i] <- xLwShellWcwc
		xInsheLLwdu[i] <- signif(xLwShellWd, 4)
		xInsheLLwdwdu[i] <- xLwShellWdwd
		xInraDu[i] <- signif(xInRad, 4) # radiation in
		xInraDwcu[i] <- xInRadwc
		xIntotaLu[i] <- signif(xIntotal, 4) # totaL exergy input
		
		## Exergy output
		xoutstorecoreu[i] <- signif(xStCore, 4) # exergy stored
		xoutstoreshelu[i] <- signif(xStShell, 4)
		xoutaIRwcu[i] <- signif(xExhaleWc, 4) # exhaled humid air
		xoutaIRwcwcu[i] <- xExhaleWcwc
		xoutaIRwdu[i] <- signif(xExhaleWd, 4)
		xoutaIRwdwdu[i] <- xExhaleWdwd
		xoutswEATwcu[i] <- signif(xSweatWc, 4) # water vapour
		xoutswEATwcwcu[i] <- xSweatWcwc
		xoutswEATwdu[i] <- signif(xSweatWd, 4) 
		xoutswEATwdwdu[i] <- xSweatWdwd
		xoutraDu[i] <- signif(xOutRad, 4) # radiation out
		xoutraDwcu[i] <- xOutRadwc
		xoutCONVu[i] <- signif(xOutConv, 4) # convection
		xoutCONVwcu[i] <- xOutConvwc
		xouttotaLu[i] <- signif(xouttotal, 4) # totaL exergy out
		
		## balance
		xconsu[i] <- signif(xConsumption, 4) # exergy consumPtion total
		tsku[i] <- signif(tsk, 4)
		tcru[i] <- signif(tcr, 4)
		wu[i] <- signif(w, 4)

	 }# End for

	#setStar <- calcet(top, tair, phia, w, im, 50, imStar); 
	results <- data.frame(
	
## Output values
		#setStar, 
		#etStar, 
		## Exergy input
		xInmetu, #metabolism 
		xInmetwcu, #metabolism warm/cold
		xInAIRwcu, 	# inhaled humid air
		xInAIRwcwcu, 
		xInAIRwdu, # 
		xInAIRwdwdu, # wet/dry
		xInLUNGwcu, # water lung
		xInLUNGwcwcu, 
		xInLUNGwdu, 
		xInLUNGwdwdu, 
		xInsheLLwcu, # water sweat
		xInsheLLwcwcu, 
		xInsheLLwdu, 
		xInsheLLwdwdu, 
		xInraDu, # radiation in
		xInraDwcu, 
		xIntotaLu, # totaL exergy input
		
		## Exergy output
		xoutstorecoreu, # exergy stored core
		xoutstoreshelu, # exergy stored shell
		xoutaIRwcu, # exhaled humid air
		xoutaIRwcwcu, 
		xoutaIRwdu, 
		xoutaIRwdwdu, 
		xoutswEATwcu, # water vapour
		xoutswEATwcwcu, 
		xoutswEATwdu, 
		xoutswEATwdwdu, 
		xoutraDu, # radiation out
		xoutraDwcu, 
		xoutCONVu, # convection
		xoutCONVwcu, 
		xouttotaLu, # totaL exergy out
		
		## balance
		xconsu, # exergy consumPtion total
		tsku, # skin temperature
		tcru, # core temperature
		wu, # skin wettedness
		stringsasFactors=FALSE
	)
	
	results
} # end of main program
###########################

