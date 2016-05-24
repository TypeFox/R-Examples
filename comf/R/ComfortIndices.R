# v1.0 done by Sophia Mueller and Marcel Schweiker

##### for maintenance #
# add variables to (1) to (3), (7) to (10)

# new function in (4) to (6)

# source("R/fctPMVPPD.r")
# source("R/fct2Node.r")
# source("R/fctdTNZ.r")
# source("R/fctATHB.r")
# source("R/fctHBxStSt.r")
# source("R/fctHBxUnSt.r")
# source("R/fctother.r")

# source("R/fcthelp.r")

calcComfInd <- function(lsCond, request="all"){

	requests = c("pmv", "ppd", "tnHumphreysNV", "tnHumphreysAC", "tAdapt15251", "dTNZ", "dTNZTa", "dTNZts", "ATHBpmv", "ATHBset", "ATHBpts", "apmv", "ptsa", "epmv", "ptse", "epCoeff", "apCoeff", "esCoeff", "asCoeff", "set", "et", "tsens", "disc", "pd", "ps", "pts", "HBxst", "pmvadj", "humidex")
	
# (1)  
  met <- as.numeric(lsCond$met)
  wme <- as.numeric(lsCond$wme)
  pb <- as.numeric(lsCond$pb)
  ltime <- as.numeric(lsCond$ltime)
  rh <- as.numeric(lsCond$rh)
  clo <- as.numeric(lsCond$clo)
  ta <- as.numeric(lsCond$ta)
  tr <- as.numeric(lsCond$tr)
  tu <- as.numeric(lsCond$tu)
  vel <- as.numeric(lsCond$vel)
  ht <- as.numeric(lsCond$ht)
  wt <- as.numeric(lsCond$wt)
  tmmo <- as.numeric(lsCond$tmmo)
  trm <- as.numeric(lsCond$trm)
  age <- as.numeric(lsCond$age)
  gender <- as.numeric(lsCond$gender) # gender: 1=female; 2=male
  tsk <- as.numeric(lsCond$tsk)
  psych <- as.numeric(lsCond$psych)
  apCoeff <- as.numeric(lsCond$apCoeff)
  epCoeff <- as.numeric(lsCond$epCoeff)
  asCoeff <- as.numeric(lsCond$asCoeff)
  esCoeff <- as.numeric(lsCond$esCoeff)
  asv <- as.numeric(lsCond$asv)
  tao <- as.numeric(lsCond$tao)
  rho <- as.numeric(lsCond$rho)
  frad <- as.numeric(lsCond$frad)
  eps <- as.numeric(lsCond$eps)
  ic <- as.numeric(lsCond$ic)
  tcrI <- as.numeric(lsCond$tcrI)
  tskI <- as.numeric(lsCond$tskI)
  basMet <- as.numeric(lsCond$basMet)
  warmUp <- as.numeric(lsCond$warmUp)
  cdil <- as.numeric(lsCond$cdil)
  sigmatr <- as.numeric(lsCond$sigmatr)

# (2)  
  l<-max(c(1, length(met), length(wme), length(pb), length(ltime), length(rh), length(clo), length(ta), length(tr), length(vel), length(ht), length(wt), length(tmmo), length(trm), length(tu), length(age), length(gender), length(vel), length(tsk), length(psych), length(apCoeff), length(epCoeff), length(asCoeff), length(esCoeff), length(asv), length(tao), length(rho), length(frad), length(eps), length(ic), length(tcrI), length(tskI), length(basMet), length(warmUp), length(cdil), length(sigmatr)))
  
  namesLsCond <- names(lsCond)

# (3)  
  paramsAll <- c("met", "wme", "pb", "ltime", "rh", "clo", "ta", "tr", "vel", "ht", "wt", "tmmo", "trm", "tu", "age", "gender", "tsk", "psych", "apCoeff", "epCoeff", "asCoeff", "esCoeff", "asv", "tao", "rho", "frad", "eps", "ic", "tcrI", "tskI", "basMet", "warmUp", "cdil", "sigmatr")
  valsAll <-c(1, 0, 760, 60, 50, 0.5, 25, 25, 0.1, 170, 70, 15, 15, 40, 21, 1, 35, (-1.4), .293, .9, (-.2), 1.3, 1.5, 5, 70, .7, .95, 1.085, 37, 36, 58.2, 60, 100, .25)
  
  # creating list of necessary input values to check for consistency and missing values
  if (length(request==1) & request[1] == "all"){
	params <- paramsAll
	vals <- valsAll
  } else {
	params <- NULL
# (4)
	for (nparam in 1:length(request)){
		if (request[nparam] == "pmv" | request[nparam] == "ppd"){
		paramsI <- c("ta", "tr", "vel", "rh", "clo", "met", "wme")
		} else if (request[nparam] == "tnHumphreysNV" | request[nparam] == "tnHumphreysAC"){
		paramsI <- c("tmmo")
		} else if (request[nparam] == "tAdapt15251"){
		paramsI <- c("trm")
		} else if (request[nparam] == "dTNZ" | request[nparam] == "dTNZTa" | request[nparam] == "dTNZTs"){
		paramsI <- c("ht", "wt", "age", "gender", "clo", "vel", "tsk", "ta", "met")
		} else if (request[nparam] == "ATHBpmv"){
		paramsI <- c("ta", "tr", "vel", "rh", "met", "wme", "psych", "trm")
		} else if (request[nparam] == "ATHBset" | request[nparam] == "ATHBpts"){
		paramsI <- c("ta", "tr", "vel", "rh", "trm", "met", "wme", "pb", "ltime", "ht", "wt", "psych")
		} else if (request[nparam] == "apmv"){
		paramsI <- c("ta", "tr", "vel", "rh", "clo", "met", "wme", "apCoeff")
		} else if (request[nparam] == "ptsa"){
		paramsI <- c("ta", "tr", "vel", "rh", "clo", "met", "wme", "pb", "ltime", "ht", "wt", "asCoeff")
		} else if (request[nparam] == "epmv"){
		paramsI <- c("ta", "tr", "vel", "rh", "clo", "met", "wme", "epCoeff", "asv")
		} else if (request[nparam] == "ptse"){
		paramsI <- c("ta", "tr", "vel", "rh", "clo", "met", "wme", "pb", "ltime", "ht", "wt", "esCoeff", "asv")			
		} else if (request[nparam] == "epCoeff" | request[nparam] == "apCoeff"){
		paramsI <- c("ta", "tr", "vel", "rh", "clo", "met", "wme", "asv")			
		} else if (request[nparam] == "esCoeff" | request[nparam] == "asCoeff"){
		paramsI <- c("ta", "tr", "vel", "rh", "clo", "met", "wme", "pb", "ltime", "ht", "wt", "asv")			
		} else if (request[nparam] == "pd"){
		paramsI <- c("ta", "tr", "vel", "rh", "clo", "met", "wme", "pb", "ltime", "ht", "wt", "tu")
		} else if (request[nparam] == "pmvadj"){
		paramsI <- c("ta", "tr", "vel", "rh", "clo", "met", "wme")
		} else if (request[nparam] == "HBxst"){
		paramsI <- c("ta", "tr", "vel", "rh", "clo", "met", "tao", "rho", "frad", "eps", "ic", "ht", "wt", "tcrI", "tskI", "basMet", "warmUp", "cdil", "sigmatr")
		} else if (request[nparam] == "humidex"){
		paramsI <- c("ta", "rh")		
		} else if (request[nparam] == "set" | request[nparam] == "et" | request[nparam] == "tsens" | request[nparam] == "disc" | request[nparam] == "ps" | request[nparam] == "pts"){
		paramsI <- c("ta", "tr", "vel", "rh", "clo", "met", "wme", "pb", "ltime", "ht", "wt")
		
		} else {
			stop(paste("error: ", request[nparam], " is not a valid request! Call listOfRequests() for all valid requests.", sep = ""))
		}
		params <- c(params, paramsI)
	}
	params <- unique(params)
	vals   <- NULL
	for (i in 1:length(params)){
		vals[i] <- valsAll[which(paramsAll == params[i])]
	}
  }
  
  for (j in 1:length(params)){ 
		# assignment of standard values for variables missing in input but required for calculation
		#if((!params[j] %in% namesLsCond) | (NA %in% get(params[j])) | (length(get(params[j]))==0)){
		if ((length(get(params[j])) > 1) & (NA %in% get(params[j]))){
			stop(paste("error: ", params[j], " is necessary for one or more of the indices required, but contains missing values!", sep=""))
		}
		if ((!params[j] %in% namesLsCond) | (NA %in% get(params[j])) | (length(get(params[j])) == 0)){
			assign(params[j], rep.int(vals[j], l))
			warning(paste("warning! ", params[j], " is necessary for one or more of the indices required, but was not given in input data. For the calculation it was set to the standard value of ", vals[j], " in all rows.", sep = ""))
		} else if (length(get(params[j])) == 1){
			assign(params[j], rep.int(get(params[j]), l))
		} else if (length(get(params[j])) != l){
			stop(paste("error: Length of ", params[j], "does not match!", sep = ""))
		}
	}
  
  #for(i in 1:l){vel[i] <- ifelse(vel[i] < 0, 0, vel[i])}
    
  if (length(request == 1) & request[1] == "all"){
      
		for (i in 1:l){
# (5)		
			giveDat <- data.frame(calcPMVPPD(ta[i], tr[i], vel[i], rh[i], clo[i], met[i], wme[i]), 
								 calc2Node(ta[i], tr[i], vel[i], rh[i], clo[i], met[i], wme[i], pb[i], ltime[i], ht[i], wt[i], tu[i], obj = "set"), 
								 calctnAuliciems(ta[i], tmmo[i]), 
								 calctnHumphreysNV(tmmo[i]), 
								 calctnHumphreysAC(tmmo[i]), 
								 calctAdapt15251(trm[i]), 
								 calcdTNZ(ht[i], wt[i], age[i], gender[i], clo[i], vel[i], tsk[i], ta[i], met[i], deltaT =.1), 
								 ATHBpmv = calcATHBpmv(trm[i], psych[i], ta[i], tr[i], vel[i], rh[i], met[i], wme[i]), 
								 ATHBset = calcATHBset(trm[i], psych[i], ta[i], tr[i], vel[i], rh[i], met[i], wme[i], pb[i], ltime[i], ht[i], wt[i]), 
								 ATHBpts = calcATHBpts(trm[i], psych[i], ta[i], tr[i], vel[i], rh[i], met[i], wme[i], pb[i], ltime[i], ht[i], wt[i]), 
								 apmv = calcaPMV(ta[i], tr[i], vel[i], rh[i], clo[i], met[i], wme[i], apCoeff[i]), 
								 epmv = calcePMV(ta[i], tr[i], vel[i], rh[i], clo[i], met[i], wme[i], epCoeff[i], asv[i]), 
								 ptsa = calcPtsa(ta[i], tr[i], vel[i], rh[i], clo[i], met[i], wme[i], pb[i], ltime[i], ht[i], wt[i], tu[i], asCoeff[i]), 
								 ptse = calcPtse(ta[i], tr[i], vel[i], rh[i], clo[i], met[i], wme[i], pb[i], ltime[i], ht[i], wt[i], tu[i], esCoeff[i], asv[i]), 
								 pmvadj = calcpmvadj(ta[i], tr[i], vel[i], rh[i], clo[i], met[i]), 
								 HBxst = calcHbExSteady(ta[i], tr[i], rh[i], vel[i], clo[i], met[i], tao[i], rho[i], frad[i], eps[i], ic[i], ht[i], wt[i], tcrI[i], tskI[i], basMet[i], warmUp[i], cdil[i], sigmatr[i])[33],
								 humidex = calcHumx(ta[i], rh[i])

								 
								 )
			if (i == 1){giveData<-giveDat}
				else{giveData<-rbind(giveData, giveDat)}
		}  
		row.names(giveData) <- 1:i 
  } else {
		for (i in 1:l){
			for (nparam in 1:length(request)){
# (6)			
				if (request[nparam] == "pmv" | request[nparam] == "ppd"){
			
					comfortData <- data.frame(calcPMVPPD(ta[i], tr[i], vel[i], rh[i], clo[i], met[i], wme[i]))
					giveDat <- with(comfortData, get(request[nparam]))
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas <- cbind(giveDatas, giveDat)}
						 
				} else if (request[nparam] == "tnHumphreysNV"){
		  
					giveDat <- calctnHumphreysNV(tmmo[i])
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas <- cbind(giveDatas, giveDat)}
						 
				} else if (request[nparam] == "tnHumphreysAC"){
		  
					giveDat <- calctnHumphreysAC(tmmo[i])
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas <- cbind(giveDatas, giveDat)}
						 
				} else if (request[nparam] == "tAdapt15251"){
		  
					giveDat <- calctAdapt15251(trm[i])
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas<-cbind(giveDatas, giveDat)}
					
				} else if (request[nparam] == "dTNZ" | request[nparam] == "dTNZTa" | request[nparam] == "dTNZTs"){
					
					comfortData <- data.frame(calcdTNZ(ht[i], wt[i], age[i], gender[i], clo[i], vel[i], tsk[i], ta[i], met[i]))		
					giveDat <- with(comfortData, get(request[nparam]))
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas<-cbind(giveDatas, giveDat)} 
					
				} else if (request[nparam] == "ATHBpmv"){
			
					giveDat <- data.frame(calcATHBpmv(trm[i], psych[i], ta[i], tr[i], vel[i], rh[i], met[i], wme[i]))
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas<-cbind(giveDatas, giveDat)}
				
				} else if (request[nparam] =="ATHBset"){
			
					giveDat <- data.frame(calcATHBset(trm[i], psych[i], ta[i], tr[i], vel[i], rh[i], met[i], wme[i], pb[i], ltime[i], ht[i], wt[i]))
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas<-cbind(giveDatas, giveDat)}
				
				} else if (request[nparam] == "ATHBpts"){
			
					giveDat <- data.frame(calcATHBpts(trm[i], psych[i], ta[i], tr[i], vel[i], rh[i], met[i], wme[i], pb[i], ltime[i], ht[i], wt[i]))
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas<-cbind(giveDatas, giveDat)}
				
				} else if (request[nparam] == "apmv"){
			
					giveDat <- data.frame(calcaPMV(ta[i], tr[i], vel[i], rh[i], clo[i], met[i], wme[i], apCoeff[i]))
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas<-cbind(giveDatas, giveDat)}

				} else if (request[nparam] == "ptsa"){
			
					giveDat <- data.frame(calcPtsa(ta[i], tr[i], vel[i], rh[i], clo[i], met[i], wme[i], pb[i], ltime[i], ht[i], wt[i], tu[i], asCoeff[i]))
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas<-cbind(giveDatas, giveDat)}
					
				} else if (request[nparam] == "epmv"){
			
					giveDat <- data.frame(calcePMV(ta[i], tr[i], vel[i], rh[i], clo[i], met[i], wme[i], epCoeff[i], asv[i]))
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas<-cbind(giveDatas, giveDat)}

				} else if (request[nparam] == "ptse"){
			
					giveDat <- data.frame(calcPtse(ta[i], tr[i], vel[i], rh[i], clo[i], met[i], wme[i], pb[i], ltime[i], ht[i], wt[i], tu[i], esCoeff[i], asv[i]))
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas<-cbind(giveDatas, giveDat)}

				} else if (request[nparam] == "pmvadj"){
			
					giveDat <- data.frame(calcpmvadj(ta[i], tr[i], vel[i], rh[i], clo[i], met[i], wme[i]))
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas<-cbind(giveDatas, giveDat)}

				} else if (request[nparam] == "epCoeff"){
					metadj <- ifelse(asv[i]>0, met[i]*(1+asv[i]*(-.067)), met[i])
					comfortData <- data.frame(calcPMVPPD(ta[i], tr[i], vel[i], rh[i], clo[i], metadj, wme[i]))
					giveDat <- with(comfortData, get("pmv"))
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas<-cbind(giveDatas, giveDat)}
					
				} else if (request[nparam] == "apCoeff"){
					comfortData <- data.frame(calcPMVPPD(ta[i], tr[i], vel[i], rh[i], clo[i], met[i], wme[i]))
					giveDat <- with(comfortData, get("pmv"))
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas<-cbind(giveDatas, giveDat)}

				} else if (request[nparam] == "esCoeff"){
					metadj <- ifelse(asv[i]>0, met[i]*(1+asv[i]*(-.067)), met[i])
					comfortData <- data.frame(calc2Node(ta[i], tr[i], vel[i], rh[i], clo[i], met[i], wme[i], pb[i], ltime[i], ht[i], wt[i], obj = "set"))
					giveDat <- with(comfortData, get("pts"))
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas<-cbind(giveDatas, giveDat)}
					
				} else if (request[nparam] == "asCoeff"){
					comfortData <- data.frame(calc2Node(ta[i], tr[i], vel[i], rh[i], clo[i], met[i], wme[i], pb[i], ltime[i], ht[i], wt[i], obj = "set"))
					giveDat <- with(comfortData, get("pts"))
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas<-cbind(giveDatas, giveDat)}
					
				} else if (request[nparam] == "pd"){
			
					comfortData <- data.frame(calc2Node(ta[i], tr[i], vel[i], rh[i], clo[i], met[i], wme[i], pb[i], ltime[i], ht[i], wt[i], tu[i], obj = "set"))
					giveDat <- with(comfortData, get(request[nparam]))
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas<-cbind(giveDatas, giveDat)}
					
				} else if (request[nparam] == "HBxst"){
					
					comfortData <- data.frame(calcHbExSteady(ta[i], tr[i], rh[i], vel[i], clo[i], met[i], tao[i], rho[i], frad[i], eps[i], ic[i], ht[i], wt[i], tcrI[i], tskI[i], basMet[i], warmUp[i], cdil[i], sigmatr[i]))		
					giveDat <- with(comfortData, get("xconss"))
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas<-cbind(giveDatas, giveDat)} 
					
				} else if (request[nparam] == "humidex"){

					giveDat <- calcHumx(ta[i], rh[i])
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas<-cbind(giveDatas, giveDat)} 
					
				} else {
			
					comfortData <- data.frame(calc2Node(ta[i], tr[i], vel[i], rh[i], clo[i], met[i], wme[i], pb[i], ltime[i], ht[i], wt[i], obj = "set"))
					giveDat <- with(comfortData, get(request[nparam]))
					if (nparam == 1){
						giveDatas<-giveDat
					} else {giveDatas<-cbind(giveDatas, giveDat)}
				}   
			}
			if (i == 1){
				giveData<-giveDatas
			} else {giveData<-rbind(giveData, giveDatas)}
	  } # end for nparam
	  
	  if (i > 1 | length(request) > 1){ #avoid warning that row.names are not unique in case of only one requested indices 
		row.names(giveData) <- 1:i 
	  }
	  giveData <- data.frame(giveData)
	  #if(length(request)>1){
	  colnames(giveData) <- request#}
	  
	} # end else
	giveData
	#comfortData
} # end fct


#create a standard data frame with all necessary variables

createCond <- function(a = TRUE){
# (7)  
  if (a == TRUE){
    ta    	<- 25   # Air temperature (degree C)
    tr    	<- 25   # mean radiant temperature (degree C)
    vel   	<- .1   # Air velocity (m/s)
    rh    	<- 50   # Relative Humidity (%)
    clo   	<- .5   # clothing (do)
    met   	<- 1.1    # metabolic rate (met)
    wme   	<- 0    # External work (met)
    tu    	<- 40   # turbulence intensity (%)
    tmmo  	<- 15   # mean monthly outdoor temperature (degree C)
    ltime 	<- 60   # Exposure time (min)
    pb    	<- 760  # Barometric pressure (torr)
    wt    	<- 70   # weight (Kg)
    ht    	<- 171  # height (cm)
	trm   	<- 15   # Running mean outdoor temperature (degree C)
	age	  	<- 21   # age (years)
	gender 	<- 1    # gender (female = 1)
	tsk   	<- 35   # mean skin temperature [degree C]
	psych 	<- (-1.4) # factor related to fixed effect on perceived control
	apCoeff <- .293 # adaptive coefficient for pmv
	epCoeff <- .9   # expectancy factor for pmv
	asCoeff <- (.2) # adaptive coefficient for set
	esCoeff <- 1.3   # expectancy factor for set
	asv     <- 1.5  # actual sensation vote (0 = neutral)
	tao		<- 5 # outdoor air temperature
	rho		<- 70   # outdoor relative humidity
	frad 	<- .7	# 0.7(for seating), 0.73(for standing) [-] 
	eps 	<- .95	# emissivity [-]
	ic 	    <- 1.085# 1.084 (average permeability), 0.4 (low permeability) 
	tcrI	<- 37	# initial values for core temp
	tskI	<- 36	# initial values for skin temperature
	basMet	<- 58.2 # basal metabolic rate
	warmUp	<- 60	# length of warm up period, i.e. number of times, loop is running for HBx calculation
	cdil	<- 100
	sigmatr <- .25

    print("A list with standard values was created. It contains standard room conditions.")
    print("For using the function 'calcComfInd', you may edit the values in that list and use this one as parameter value or create one on your own.")

# (8)    
    list(ta = ta, tr = tr, vel = vel, rh = rh, clo = clo, met = met, wme = wme, tu = tu, tmmo = tmmo, ltime = ltime, pb = pb, wt = wt, ht = ht, trm = trm, age = age, gender = gender, tsk = tsk, psych = psych, apCoeff = apCoeff, epCoeff = epCoeff, asCoeff = asCoeff, esCoeff = esCoeff, asv = asv, tao = tao, rho = rho, frad = frad, eps = eps, ic = ic, tcrI = tcrI, tskI = tskI, basMet = basMet, warmUp = warmUp, cdil = cdil, sigmatr = sigmatr)
    
  } else {
# (9)
    ta    	<- NA # Air temperature (degree C)
    tr    	<- NA # mean radiant temperature (degree C)
    vel   	<- NA # Air velocity (m/s)
    rh    	<- NA # Relative Humidity (%)
    clo   	<- NA # clothing (do)
    met   	<- NA # metabolic rate (met)
    wme   	<- NA # External work (met)
    tu    	<- NA # turbulence intensity (%)
    tmmo  	<- NA # mean monthly outdoor temperature (degree C)
    ltime 	<- NA # Exposure time (min)
    pb    	<- NA # Barometric pressure (torr)
    wt    	<- NA # weight (Kg)
    ht    	<- NA # height (cm)
	trm   	<- NA # Running mean outdoor temperature (degree C)
	age	  	<- NA # age (years)
	gender 	<- NA # gender (female = 1)
	tsk   	<- NA # mean skin temperature [degree C]
	psych 	<- NA # factor related to fixed effect on perceived control
	apCoeff  <- NA # adaptive coefficient
	epCoeff  <- NA # expectancy factor
	asCoeff  <- NA # adaptive coefficient for set
	esCoeff  <- NA  # expectancy factor for set
	asv     <- NA # actual sensation vote (0 = neutral)
	tao		<- NA # outdoor air temperature
	rho		<- NA   # outdoor relative humidity
	frad 	<- NA	# 0.7(for seating), 0.73(for standing) [-] 
	eps 	<- NA	# emissivity [-]
	ic 	    <- NA   # 1.084 (average permeability), 0.4 (low permeability) 
	tcrI	<- NA	# initial values for core temp
	tskI	<- NA	# initial values for skin temperature
	basMet	<- NA   # basal metabolic rate
	warmUp	<- NA	# length of warm up period, i.e. number of times, loop is running for HBx calculation
	cdil	<- NA
	sigmatr <- NA

	print("An empty list was created. You need to edit it for using it in the function 'calcComfInd'.")
 
#(10)
    list(ta = ta, tr = tr, vel = vel, rh = rh, clo = clo, met = met, wme = wme, tu = tu, tmmo = tmmo, ltime = ltime, pb = pb, wt = wt, ht = ht, trm = trm, age = age, gender = gender, tsk = tsk, psych = psych, apCoeff = apCoeff, epCoeff = epCoeff, asCoeff = asCoeff, esCoeff = esCoeff, asv = asv, tao = tao, rho = rho, frad = frad, eps = eps, ic = ic, tcrI = tcrI, tskI = tskI, basMet = basMet, warmUp = warmUp, cdil = cdil, sigmatr = sigmatr)
    
 }
}

listOfRequests <- function(){
	print(c("pmv", "ppd", "tnHumphreysNV", "tnHumphreysAC", "tAdapt15251", "dTNZ", "dTNZTa", "dTNZts", "ATHBpmv", "ATHBset", "ATHBpts", "apmv", "ptsa", "epmv", "ptse", "epCoeff", "apCoeff", "esCoeff", "asCoeff", "set", "et", "tsens", "disc", "pd", "ps", "pts", "HBxst", "pmvadj", "humidex"))
}