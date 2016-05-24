# Functions return:
#
# dTNZ - 
# dTNZTa - 
# dTNZTs - 
#
# File contains 1 function:
#   - calcdTNZ(182, 82, 21, 1, .5, .1, 36, 20)
#       returns dTNZ, dTNZTa, dTNZTs
#
# v1.0 done by Marcel Schweiker in cooperation with B. Kingma


# Function: dTNZ ################
###########################################
		
calcdTNZ <- function(ht, wt, age, gender, clo, vel, tsk, ta, met, deltaT =.1){		

	# calculation of surface area according to duBois
		#sa <- (wt ^ .425*ht ^ .725)*.007184
	# calculation of surface area according to mosteller
	sa <- (wt ^ .5 * ht ^ .5) * 1 / 60

	# calculation of basal metabolic rate according to Harris-Benedict equations revised by Roza and Shizgal in 1984
	if (gender == 1){ # female
		basMet <- 447.593 + (9.247 * wt) + (3.098 * ht) - (4.330 * age)
	} else {
		basMet <- 88.362 + (13.397 * wt) + (4.799 * ht) - (5.677 * age)
	}
	basMet <- 4164 / 86400 * basMet
	basMet <- basMet * 5 / 4 # adjust for basMet = .8 met

	# definition of resolution of combinations
	deltaT  <- .1
	mult    <- 1 / deltaT
	seqi    <- seq(deltaT, 20, deltaT)
	seqj    <- seq(deltaT, 26, deltaT)
	offseqi <- 19
	offseqj <- 9

	alpha  <- 0.08 # []
	w      <- 0.06 # skin wettedness fraction []
	gammac <- 0.00750061683 # [mmHg/pa]
	lambda <- 2.2 # [degree C/mmHg] Lewis relation
	phi    <- 0.50 # [%]
	Tcmin  <- 36 # [degree C]
	Tcmax  <- 38 # [degree C]

	cloV <- clo
	velV  <- vel
	velVStill <- 100 * 0.05
	tskObs <- tsk
	taObs  <- ta

	A <- sa
	dTNZHori <- wert <- dTNZVert <- dTNZ <- dTNZTs <- dTNZTa <- NA
	# get TNZ for each timestept of data
	#for ( k in 1:nrow(dfsubj)){

	icl <- 0.155 * cloV # [m2 degree C/W]
	va <- velV # [m/s]
	# convert va to cm/s
	va <- va * 100 # [cm/s]

	mmin <- basMet * (met-.1) # [W]
	mmax <- basMet * (met+.1) # [W]

	IBodymax <- 0.112 # [m2 degree C/W] see Veicsteinas et al.
	IBodymin <- 0.031 # [m2 degree C/W]
	Tsmin <- Tcmin - (1 - alpha) * mmax * IBodymax / A # [degree C] 
	Tsmax <- Tcmax - (1 - alpha) * mmin * IBodymin / A # [degree C] 

	# define empty matrix
	Tc <- IBodyx <- tax <- Tsx <- Iairx <- Q <- data.frame(matrix(ncol = length(seqi), nrow = length(seqj)))
	
	ta <- matrix(rep((seqj + offseqj), length(seqi)), length(seqj), length(seqi))
	Ts <- matrix(rep((seqi + offseqi), each = length(seqj)), length(seqj), length(seqi))

	IBody <- IBodymax + ((IBodymax - IBodymin) / (Tsmin - Tsmax)) * (Ts - Tsmin)
	IBody <- ifelse(IBody > IBodymax, IBodymax, IBody)
	IBody <- ifelse(IBody < IBodymin, IBodymin, IBody)

	Iair <- 1 / ((0.19 * sqrt(va) * (298 / (ta + 273.15))) + (0.61 * ((ta + 273.15) / 298) ^ 3)) 
	#a1 <- 0.8*(0.19*sqrt(velVStill)*(298/(ta+273.15))) + 0.2*(0.19*sqrt(va)*(298/(ta+273.15)))
	#a2 <- (0.61*((ta+273.15)/298) ^ 3)
	#Iair <- 1/(a1+a2) # weighted average due to body not completely exposed to air velocity

	Iair <- Iair * 0.155 # to adjust from clo unit to m2K/w

	hconv <- 1 / Iair - (0.61 * ((ta + 273.15) / 298) ^ 3) / 0.155 # [W/m2 degree C]
	#hconv <- 0.19*sqrt(va)*(298/(ta+273.15))

	ps <- gammac * 100 * exp(18.965-4030 / (Ts + 235))
	pair <- gammac * phi * 100 * exp(18.965 - 4030/(ta + 235))

	fpcl <- 1 / (1 + hconv * icl)

	Qe  <- A * w * lambda * hconv * (ps - pair) * fpcl
	Qrc <- (A / (icl + Iair)) * (Ts - ta) 

	# calculating and storing values for Tc
	Tc <- Ts + (IBody / (1 - alpha)) * ((Ts - ta) / (icl + Iair) + (Qe / A))
	# calculation and storing values for Q: metabolic heat production
	Q <- Qrc + Qe
	## storing IBody values in matrix
	#IBodyx[jmult, imult] <- IBody
	## storing Ambient temperature values in matrix
	tax <- ta
	## storing skin temperature values in matrix
	Tsx <- Ts
	## storing Iair values in matrix
	#Iairx[jmult, imult] <- Iair
	 
	colnames(Tc) <- colnames(Q) <- seqi + offseqi # 
	rownames(Tc) <- rownames(Q) <- seqj + offseqj #

	#setting all to NA, where Tc is unrealistic 36 < Tc < 38
	Q[Tc > Tcmax] <- NA
	Q[Tc < Tcmin] <- NA

	Tc[Tc > Tcmax] <- NA
	Tc[Tc < Tcmin] <- NA

	# defintion of vectors used for image()
	taGr <- seqj + offseqj 
	TsGr <- seqi + offseqi 

	Qtnz <- Q
	Qtnz[Q < (1 - alpha) * mmin] <- NA
	Qtnz[Q > (1 - alpha) * mmax] <- NA

	#########################
	# from here get dTNZ (vertical + horizontal)
	#########################

	# transfer values in Qtnz to values of columns (tsk)
	QtnzTs <- Qtnz

	listTs <- NA
	g <- 1

	for (i in 1:ncol(Qtnz)){
		for (j in 1:nrow(Qtnz)){
			if (!is.na(Qtnz[j, i])){
				# QtnzTs[j, i] <- as.numeric(colnames(Qtnz)[i])
				 listTs[g] <- as.numeric(colnames(Qtnz)[i])
				 g <- g + 1
			 }
		 }
	 }
	#QtnzTs
	tskCentroid <- median(listTs)

	# transfer values in Qtnz to values of columns (tsk)
	Qtnztop <- Qtnz

	listtop <- NA
	g <- 1

	for (i in 1:ncol(Qtnz)){
		for (j in 1:nrow(Qtnz)){
			if (!is.na(Qtnz[j, i])){
				# Qtnztop[j, i] <- as.numeric(row.names(Qtnz)[j])
				 listtop[g] <- as.numeric(row.names(Qtnz)[j])
				 g <- g + 1
			 }
		 }
	 }
	#Qtnztop
	topCentroid <- median(listtop)

	#%calculate distance to TNZ (dTNZ)
	dTNZ <- round(sqrt((taObs-topCentroid) ^ 2 + (tskObs-tskCentroid) ^ 2 ), 2); # abs value of distance
	dTNZTs <- round(tskObs-tskCentroid, 2) # rel value of distance assuming tskin to be dominant for sensation
	dTNZTa <- round(taObs-topCentroid, 2) # rel value of distance assuming tambient to be dominant for sensation
	tskCentroid <- round(tskCentroid, 2)
	topCentroid <- round(topCentroid, 2)
	
	data.frame(dTNZ, dTNZTs, dTNZTa, tskCentroid, topCentroid)

	}
		
		
