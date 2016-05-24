lhc_generate_lhc_sample <-
function(FILEPATH,PARAMETERS,NUMSAMPLES,PMIN,PMAX,ALGORITHM)
{
	if(requireNamespace("lhs",quietly=TRUE))
	{
		if(file.exists(FILEPATH))
		{
			NUMPARAMS<-length(PARAMETERS)

			ALGORITHM<-tolower(ALGORITHM)

			# PERFORM THE SAMPLING - JUDGING ON USERS CHOICE OF ALGORITHM
			if(ALGORITHM=="optimum")
			{
				# MAY TAKE A WHILE FOR A LARGE NUMBER OF SAMPLES (THIS BEING TWO DAYS WHERE NUMSAMPLES=500)
				design<-lhs::optimumLHS(NUMSAMPLES,NUMPARAMS,maxSweeps=2,eps=.1)
			}
			else
			{
				design<-lhs::randomLHS(NUMSAMPLES,NUMPARAMS)
			}


			# NOW LOOK AT THE VALUE CHOSEN FOR EACH SAMPLE, AS THESE WILL CURRENTLY BE BETWEEN 0 AND 1 
			for(k in 1:NUMSAMPLES)
			{
				# NOW LOOK AT EACH PARAMETER IN TURN
				# THE LHC WILL HAVE GIVEN EACH A VALUE BETWEEN 0 AND 1
				# NOW USE THE MAX AND MIN VALUES FOR EACH PARAMETER TO GIVE IT A PROPER VALUE
				for(l in 1:NUMPARAMS)
				{
					# GET THE MAX AND MIN VALUES FOR THIS PARAMETER FROM THE ARRAY
					lhc_max<-PMAX[l]
					lhc_min<-PMIN[l]
		
					# NOW CALCULATE THE VALUE TO USE FOR THIS PARAMETER
					value<-(design[k,l] * (lhc_max - lhc_min) ) + lhc_min
		
					# NOW REPLACE THE VALUE IN THE TABLE (BETWEEN 0 AND 1) WITH THE PARAMETER VALUE
					design[k,l]<-value          
				}
			}                     

			# LABEL THE RESULTS
			colnames(design)<-c(PARAMETERS)

			# OUTPUT THE LHC DESIGN AS A CSV FILE
			write.csv(design,paste(FILEPATH,"/LHC_Parameters_for_Runs.csv",sep=""),row.names=FALSE,quote = FALSE)

			print(paste("Parameter Set Generated and Output to ",FILEPATH,"/LHC_Parameters_for_Runs.csv",sep=""))
		}
		else
		{
			print("The directory specified in FILEPATH does not exist. No parameter samples generated")
		}
	}
	else
	{
		print("The lhc_generate_lhc_sample function requires the lhs package to be installed")
	}
}

