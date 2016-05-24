Lumped_VSA_model <-
function(
	dateSeries, # Date series (format should be: "2003-02-36")
	P, # Rain + Snow melt (mm)
	Tmax, # Max daily temperature (C)
	Tmin, # Max daily temperature (C)
	Depth = NULL, # Average soil depth in watershed (mm) [don't need if AWC and SAT entered directly]
	SATper = NULL, # Porosity (fraction)
	AWCper = NULL, # Available Water Capacity (AWC) (fraction)
	percentImpervious = 0, # Percent of the watershed that is impervious
	no_wet_class = 10, # The number of wetness classes for saturated area designation
	Tp = 5, # Time to peak (hours)
	latitudeDegrees = 42.38, 
	albedo = 0.23, # Average albedo
	StartCond = "avg", # Watershed conditions before first day of run ("wet", "dry", "avg")
	PETin = NULL,# User has the option to enter PET values (mm/day)
	AWC = Depth*AWCper, # AWC as a depth (mm)
	SAT = Depth*SATper,# porosity as a depth (mm)
	SW1 = NULL, # Soil water on the first day (depth, mm)
	BF1 = 1, #  mm/day can use nearby watershed baseflow, ConvertFlowUnits(cfs,WA=W_Area)
	PETcap = 5, # Does not let PET get larger than this cut off (mm)
	rec_coef = 0.1, # based on a study by Weiler in NY state
	Se_min = 78, # mm
	C1 = 3.1, # Coefficient relating soil water to Curve Number S
	Ia_coef = 0.05, # range ~ 0.05 - 0.2
	PreviousOutput = NULL, # Allows us to take previous model output to initiate the model
	runoff_breakdown = RunoffBreakdown(Tp, HrPrcDelay = (Tp/2-4)) 
	  # The proportion of runoff that reaches the outlet on a given day after the storm event. 
	  # Calculated from Time to peak, which is related to time of concentration
){
	# Parameters
	latitude<-latitudeDegrees*pi/180 ## latitude in radians

	## Effective soil water storage coefficients, eswsc = Se*(2..... see Schneiderman et al 2007)
	eswsc <- vector(mode="numeric", length=no_wet_class)
	eqn16 <- function(x){(2*(sqrt(1-(1/no_wet_class)*(x-1))-sqrt(1-(1/no_wet_class)*x))/(1/no_wet_class))-1}
	x <- seq(1,no_wet_class, by=1)
	eswsc <- eqn16(x)


	### This is to allow us to use previous model output as an input - useful for calculating just
	### a few more days of modeled flow without having to run years of data ..
	if (!is.null(PreviousOutput)){# We will extend the input variables to include overlap with previous model input
		nr<- nrow(PreviousOutput)
		dateSeries<- c(PreviousOutput$Date[(nr-3):nr],dateSeries)
		P<- c(PreviousOutput$rain_snowmelt[(nr-3):nr], P)
		Tmax<- c(PreviousOutput$Tmax[(nr-3):nr], Tmax)
		Tmin<- c(PreviousOutput$Tmin[(nr-3):nr], Tmin)
		OutStart<- 5# The output that we report will not include previous output
	} else OutStart<- 1
	################################################
	Tav<-(Tmax+Tmin)/2
	day<-strptime(dateSeries,"%Y-%m-%d")$yday+1
	month<-strptime(dateSeries,"%Y-%m-%d")$mon+1

	## Potential Evapotranspiration 
	PET<-PET_fromTemp(Jday=day,Tmax_C=Tmax,Tmin_C=Tmin,AvgT=Tav,albedo=albedo,lat_radians=latitude)*1000## mm (Priestley-Taylor)
	PET[which(PET>PETcap)]<-PETcap#Sets a cap on PET estimates (usually ~ 5mm)
	
	ETo <- PET
	ETo[which(day<=166)]<-(0.1+0.02*(day[which(day<=166)]-121))*PET[which(day<=166)]## linear increase from May 1- June 15
	ETo[which(month<5)]<-0.1*PET[which(month<5)]									## until May 1, ET is only 10% of PET
	ETo[which(day>=274)]<-(1-0.02*(day[which(day>=274)]-274))*PET[which(day>=274)]  # 
	ETo[which(day>319)]<-0.1*PET[which(day>319)]									## After Nov 15, ET is only 10% of PET

	if (!is.null(PETin)) ETo <- PETin  ## Allows users to insert PET themselves
	ET <- ETo    #  Initializing modeled actual ET, equal to ETo when deltaP > 0, but less than ETo on drier days
	# Values to be defined in a loop

	# Initializing vectors for Water Budget Loop, and Runoff estimation

	SoilWater<-vector(length=length(P))##(mm)
	excess<-vector(length=length(P))
	TM_S<-vector(length=length(P))##This is the daily time-step T-M storage, used to calc baseflow
	totQ<-vector(length=length(P))
	Se<-vector(length=length(P))
	Q2<-vector(length=length(P))
	baseflow<-vector(length=length(P))
	MaxWetClass<-vector(length=length(P))
	impervRunoff<-vector(length=length(P))
	OverlandFlow<-vector(length=length(totQ))
	ShallowInterflow<-vector(length=length(totQ))

	deltaP<-P - ETo## (mm) neg values indicate net evap
	impervIa<-Ia_coef * 5  ##Se = 5 mm for impervious areas (CN=98)
	impervPe<-deltaP - impervIa##  Effective precipitation on impervious areas
	impervPe[which(impervPe < 0)]<-0
	Pe<-vector(length=length(P))# Effective Precipitation (non-impervious areas)
	sigmaS<-matrix(nrow=length(P), ncol=no_wet_class)  # Storage in each wetness class
	runoff_by_class<-matrix(nrow=length(P), ncol=no_wet_class)


	# Setting up initial values
	TM_S[1]<-BF1/rec_coef
	if(OutStart==5){# If previous model data is used, initialize with that
		SoilWater[1:4]<-PreviousOutput$SoilWater[(nr-3):nr]
		excess[1:4]<-PreviousOutput$excess[(nr-3):nr]
		TM_S[1:4]<-PreviousOutput$baseflow[(nr-3):nr]/rec_coef#baseflow_obs[1]/rec_coef
		totQ[1:4]<-PreviousOutput$totQ[(nr-3):nr]
		OverlandFlow[1:4]<-PreviousOutput$OverlandFlow[(nr-3):nr]
		ShallowInterflow[1:4]<-PreviousOutput$ShallowInterflow[(nr-3):nr]
		Se[1:4]<-PreviousOutput$Se[(nr-3):nr]
		sigmaS[1,]<-eswsc * Se[(nr-3)]
		sigmaS[2,]<-eswsc * Se[nr-2]
		sigmaS[3,]<-eswsc * Se[nr-1]
		sigmaS[4,]<-eswsc * Se[nr]
	} else {# Otherwise we will initialize with average values based on initial conditions ("wet", "dry", "avg")
		if (is.null(SW1)){
			SoilWater[1] <- ifelse(StartCond == "wet", AWC, ifelse(StartCond == "dry", 0.2*AWC, 0.65*AWC))
		} else SoilWater[1] <- SW1# If user knows soil water on first day, use that
		excess[1] <- 0# Assume no runoff from the day before
		Se[1]<- Se_min+C1*(AWC-SoilWater[1])
		sigmaS[1,] <- eswsc*Se[1]
		totQ[1] <- Pe[1]*Pe[1]/(Pe[1]+Se[1])
		impervRunoff[1] <- impervPe[1]
	}
	Ia<-Ia_coef*Se  # Initial abstraction = depth of precip before runoff starts

	## Thornthwaite-Mather Function and Runoff Generation 
	for(i in (2):length(deltaP)){  
		# First calculate runoff
		if (deltaP[i]-Ia[i-1]>0){
			Pe[i] <- deltaP[i]-Ia[i-1]  
		} else Pe[i]<-0

		totQ[i] <- Pe[i]*Pe[i]/(Pe[i]+Se[i-1]) ## Effective storage is from previous day

		# Water Balance
		if(deltaP[i]<=0){  #DRYING CONDITION
		SoilWater[i]<- SoilWater[i-1]*exp(deltaP[i]/AWC)
		ET[i] <- SoilWater[i-1]*(1-exp(deltaP[i]/AWC))  ## amount that gets evaporated (mm)
		} else if (deltaP[i]+SoilWater[i-1] <= AWC){  # WETTING BELOW FIELD CAPACITY
		SoilWater[i]<- deltaP[i]+SoilWater[i-1]
		} else { # WETTING ABOVE FIELD CAPACITY
		SoilWater[i] <- AWC  ## So overall SW cannot exceed AWC
		excess[i] <- deltaP[i]+SoilWater[i-1]-AWC
		}
		Se[i] <- Se_min+C1*(AWC-SoilWater[i-1])

		sigmaS[i,]<- eswsc*Se[i]  ## Amount of storage available in each wetness class [mm]
		Ia[i]<-Ia_coef*Se[i]

		if ((excess[i])>=totQ[i]){# Ensure mass-balance, since curve number is empirical
		TM_S[i]<-TM_S[i-1]*(1-rec_coef)+excess[i]-totQ[i]  
		} else {
		TM_S[i]<-TM_S[i-1]*(1-rec_coef)
		SoilWater[i]<-SoilWater[i]-(totQ[i]-excess[i])
		}

		baseflow[i]<-rec_coef*TM_S[i]
		impervRunoff[i]<-impervPe[i] 

	}
	### End of Water Balance Loop ############################################################################

	# Use the coefficients generated from time of concentration
	runoff_breakdown[which(runoff_breakdown < 0.01)] <- 0
	if (OutStart==1){  # Initialize these values if not already taken from previous output
	OverlandFlow[1:4]<-totQ[1:4]*runoff_breakdown[1]
	ShallowInterflow[1]<-0# Assuming no runoff in previous days before this model run
	ShallowInterflow[2]<-totQ[1]*runoff_breakdown[2]
	ShallowInterflow[3]<-totQ[1]*runoff_breakdown[3] + totQ[2]*runoff_breakdown[2]
	ShallowInterflow[4]<-totQ[1]*runoff_breakdown[4] + totQ[2]*runoff_breakdown[3] + totQ[3]*runoff_breakdown[2]
	}

	for (i in 5:length(totQ)){

	OverlandFlow[i]<- totQ[i]*runoff_breakdown[1]
	ShallowInterflow[i]<- totQ[i-4]*runoff_breakdown[5]+totQ[i-3]*runoff_breakdown[4]+totQ[i-2]*runoff_breakdown[3]+totQ[i-1]*runoff_breakdown[2]

	#Determine the number of wetness classes contributing to flow
	if (OverlandFlow[i]+ShallowInterflow[i]>0){
	MaxRunoffOnlyInEachClass<-vector(length=no_wet_class)
	for(j in 1:(no_wet_class-1)){
	MaxRunoffOnlyInEachClass[j]<-sigmaS[i,j+1]-sigmaS[i,j]
	}
	MaxRunoffOnlyInEachClass[no_wet_class]<-100 ## Last wetness class is assigned a high value - will not generate runoff.
	intermediateSum<-0
	j <- 1  # Now cycle through and determine MaxWetClass contributing to flow
	while (j < (no_wet_class+1) ){
	if (((MaxRunoffOnlyInEachClass[j]*j+intermediateSum)/no_wet_class)>=OverlandFlow[i]+ShallowInterflow[i]){
	MaxWetClass[i]<- j   ## so MaxWetClass represents last wetness class that produced runoff
	runoff_by_class[i,j]<- ((OverlandFlow[i]+ShallowInterflow[i])*no_wet_class-intermediateSum)/j
	for (k in (j-1):1){
	runoff_by_class[i,k] <- sigmaS[i,j]-sigmaS[i,k]+runoff_by_class[i,j]
	}
	j<-no_wet_class+1  ## ends the while loop
	} else {
	intermediateSum<-sum(sigmaS[i,j+1]-sigmaS[i,1:j])
	j<-j+1
	}
	}
	} else MaxWetClass[i] <- 0 # If no runoff, then no saturated areas.
	}
	MaxWetPer <- MaxWetClass/no_wet_class
	
	modeled_flow<-(OverlandFlow+ShallowInterflow)*(1-0.01*percentImpervious)+impervRunoff*(0.01*percentImpervious)+baseflow#totQ+baseflow
	quickflow_combined<-(OverlandFlow+ShallowInterflow)*(1-0.01*percentImpervious)+impervRunoff*(0.01*percentImpervious)
	rain_snowmelt<-P
	Streamflow<-data.frame(Date=as.Date(as.character(dateSeries)), rain_snowmelt, modeled_flow, baseflow, OverlandFlow, ShallowInterflow, totQ, deltaP, Pe, Se, Ia, SoilWater, quickflow_combined, impervRunoff, excess, Tmax, Tmin, ET,MaxWetClass,MaxWetPer)[OutStart:length(P),] 

	return(Streamflow)
}
