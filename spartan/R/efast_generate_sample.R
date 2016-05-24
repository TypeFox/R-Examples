efast_generate_sample <-
function(FILEPATH,NUMCURVES,NUMSAMPLES,PARAMETERS,PMIN,PMAX)
{
	if(file.exists(FILEPATH))
	{
		NUMPARAMS<-length(PARAMETERS)
		wantedN <- NUMSAMPLES*NUMPARAMS*NUMCURVES # wanted no. of sample points

		# OUTPUT
		# SI[] : first order sensitivity indices
		# STI[] : total effect sensitivity indices
		# Other used variables/constants:
		# OM[] : vector of k frequencies
		# OMi : frequency for the group of interest
		# OMCI[] : set of freq. used for the compl. group
		# X[] : parameter combination rank matrix
		# AC[],BC[]: fourier coefficients
		# FI[] : random phase shift
		# V : total output variance (for each curve)
		# VI : partial var. of par. i (for each curve)
		# VCI : part. var. of the compl. set of par...
		# AV : total variance in the time domain
		# AVI : partial variance of par. i
		# AVCI : part. var. of the compl. set of par.
		# Y[] : model output
	
		MI <- 4  # maximum number of fourier coefficients
		# that may be retained in calculating the partial
		# variances without interferences between the
		# assigned frequencies
	
		# Computation of the frequency for the group
		# of interest OMi and the # of sample points NUMSAMPLES (here N=NUMSAMPLES)
		OMi <- floor(((wantedN/NUMCURVES)-1)/(2*MI)/NUMPARAMS)
		NUMSAMPLES <- 2*MI*OMi+1
		if(NUMSAMPLES*NUMCURVES < 65)
		{
		    print("Error: sample size must be >= 65 per factor")
		}

		PARAMETERVALS<-array(0,dim=c(NUMSAMPLES,NUMPARAMS,NUMPARAMS,NUMCURVES))

		for(PARAMNUM in 1:NUMPARAMS) #loop through parameters
		{
			# Algorithm for selecting the set of frequencies.
			# OMci(i), i=1:k-1, contains the set of frequencies
			# to be used by the complementary group.
	
			OMci <- efast_setfreq(NUMPARAMS,OMi/2/MI,PARAMNUM)
			OM<-array(0,dim=c(1,NUMPARAMS,1))
	
		    	# Loop over the NUMCURVES search curves.
			for(CURVENUM in 1:NUMCURVES)
		    	{
		        	# Setting the vector of frequencies OM
		        	# for the k parameters
		        	cj <- 1
				for(j in 1:NUMPARAMS)
				{
		                	if(j==PARAMNUM)
					{
		                		# For the parameter (factor) of interest
						# RECODE WHEN WORKED OUT OM ARRAY                		
						OM[PARAMNUM] = OMi;
					}else
					{
		                		# For the complementary group.
						# RECODE WHEN WORKED OUT ARRAY
		                		OM[j] = OMci[cj]
		                		cj <- cj+1
					}
				}
	          
		       		# Setting the relation between the scalar
		        	# variable S and the coordinates
		        	# {X(1),X(2),...X(k)} of each sample point.
				FI<-array(runif(NUMPARAMS,min=0,max=1),dim=c(NUMPARAMS,1,1))
				FI<-FI*2*pi
		        
				S_VEC <- pi*(2*(1:NUMSAMPLES)-NUMSAMPLES-1)/NUMSAMPLES
		        	OM_VEC <- OM[1:NUMPARAMS]
		
				FI_MAT <- array(0,dim=c(NUMPARAMS,NUMSAMPLES,1))
			
				for(i in 1:NUMSAMPLES)
				{
					FI_MAT[,i,1]<-FI
				}
		
				# FORMULA IN ORIGINAL MATLAB CODE:
		        	#ANGLE = OM_VEC'*S_VEC+FI_MAT;
				# CONVERSION TO R:
				omVecSVec<-array(OM_VEC%*%t(S_VEC),dim=c(NUMPARAMS,NUMSAMPLES,1))
				ANGLE<-omVecSVec+FI_MAT
	
				# TRANSPOSE ARRAY
				ANGLET<-array(0,dim=c(NUMSAMPLES,NUMPARAMS,1))
				for(i in 1:NUMSAMPLES)
				{
					ANGLET[i,,1]<-ANGLE[,i,1]
				}
			
				# NOW CALCULATE THE PARAMETER VALUES - THESE ARE STORED IN A MULTIDIMENSIONAL ARRAY, AS EACH CURVE HAS SEVEN SETS OF PARAMETER VALUES
		        	PARAMETERVALS[,,PARAMNUM,CURVENUM] <- 0.5+asin(sin(ANGLET))/pi
	       		
				# AS THESE VALUES WILL CURRENTLY BE BETWEEN 0 AND 1, TRANSFORM THE DISTRIBUTION TO GIVE TRUE PARAMETER VALUES
				PARAMETERVALS[,,PARAMNUM,CURVENUM] <- efast_parameterdist(PARAMETERVALS[,,PARAMNUM,CURVENUM],PMAX,PMIN,NUMSAMPLES,length(PARAMETERS))
			}
		}

		# NOW OUTPUT THE RESULTS - SPLIT BY CURVE FILE
		# SO, WILL HAVE ONE FILE FOR EACH PARAMETER OF INTEREST, FOR EACH CURVE
		for(CURVENUM in 1:NUMCURVES)
		{
			for(PARAMNUM in 1:NUMPARAMS)
			{
				parameterFile = paste(FILEPATH,"/Curve",CURVENUM,"_",PARAMETERS[PARAMNUM],".csv",sep="")
				outputParams <- PARAMETERVALS[,,PARAMNUM,CURVENUM]
				colnames(outputParams)<-c(PARAMETERS)
				
				write.csv(outputParams,parameterFile,quote = FALSE,row.names=FALSE)

				print(paste("Parameter Set for ",CURVENUM," Generated and Output to ",FILEPATH,"/Curve",CURVENUM,"_",PARAMETERS[PARAMNUM],".csv",sep=""))
			}
		}
		

	}
	else
	{
		print("The directory specified in FILEPATH does not exist. No parameter samples generated")
	}
}

