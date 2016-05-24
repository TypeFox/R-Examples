efast_sd <-
function(RESULTSARRAY,OMi,MI,OUTMEASURES,NUMPARAMS,NUMCURVES)
{
	# RESULTSMAT FORMAT - FIRST: PARAM SET NUMBER,2ND: TIMEPOINT (NOT USED BUT IS IN EXAMPLE),3RD: RESULT VALUE ARRAY,4TH:PARAMETER,5TH: CURVE
	#[a b c NUMPARAMS NUMCURVES]=size(RESULTSMAT);
	#if nargin<5
	#    display(['ERROR = Choose one or more outputs from var 1 and variable ',num2str(c),' of the model'])
	#    error('eFAST: the output components for the sensitivity is missing. Not enough input arguments.')

	Si<-array(0,dim=c(NUMPARAMS,1,OUTMEASURES))
	STi<-array(0,dim=c(NUMPARAMS,1,OUTMEASURES))
	rangeSi<-array(0,dim=c(NUMPARAMS,NUMCURVES,OUTMEASURES))
	rangeSTi<-array(0,dim=c(NUMPARAMS,NUMCURVES,OUTMEASURES))
	Vci<-array(0,dim=c(1,NUMCURVES,1))
	Vi<-array(0,dim=c(1,NUMCURVES,1))
	V<-array(0,dim=c(1,NUMCURVES,1))

	for(MEASURE in 1:OUTMEASURES)
	{
		for(PARAMNUM in 1:NUMPARAMS) #loop through parameters
		{
            		# Initialize AV,AVi,AVci to zero.
			# THOUGH THESE SEEM TO BE HIGHLIGHTED OUT OF THE ORIGINAL
            		AV<-0;
            		AVi<-0;
            		AVci<-0;

			for(CURVENUM in 1:NUMCURVES)  
			{
				# GET THE RESULTS FOR THIS CURVE, FOR THIS PARAMETER, FOR THIS MEASURE
				# THEN SUBTRACT THE MEAN OF THE COLUMN FROM THE OUTPUT VALUE
				#print(RESULTSARRAY[,(((PARAMNUM*OUTMEASURES)-OUTMEASURES)+MEASURE),CURVENUM])
				
				MEASURE_RESULTS_FOR_PARAM<-na.omit(RESULTSARRAY[,(((PARAMNUM*OUTMEASURES)-OUTMEASURES)+MEASURE),CURVENUM])
				MEASURE_RESULTS_FOR_PARAM<-MEASURE_RESULTS_FOR_PARAM-t(mean(MEASURE_RESULTS_FOR_PARAM))

				#RESULTSARRAY[,(((PARAMNUM*OUTMEASURES)-OUTMEASURES)+MEASURE),CURVENUM]<-
				#	RESULTSARRAY[,(((PARAMNUM*OUTMEASURES)-OUTMEASURES)+MEASURE),CURVENUM] - 
				#	t(mean(RESULTSARRAY[,(((PARAMNUM*OUTMEASURES)-OUTMEASURES)+MEASURE),CURVENUM]))
                
				# Fourier coeff. at [1:OMi/2].
				# GET THE NUMBER OF SAMPLES FOR THIS OUTPUT
				N<-length(MEASURE_RESULTS_FOR_PARAM)
				#N<-length(RESULTSARRAY[,(((PARAMNUM*OUTMEASURES)-OUTMEASURES)+MEASURE),CURVENUM])
				# NQ GOES JUST BELOW THE MIDPOINT
                		NQ<-(N-1)/2
				# NO GOES JUST ABOVE THE MIDPOINT
                		N0<-NQ+1
                		
				Y_VECP <- MEASURE_RESULTS_FOR_PARAM[N0+(1:NQ)] + MEASURE_RESULTS_FOR_PARAM[N0-(1:NQ)]

				Y_VECM <- MEASURE_RESULTS_FOR_PARAM[N0+(1:NQ)] - MEASURE_RESULTS_FOR_PARAM[N0-(1:NQ)]

				#Y_VECP = RESULTSARRAY[N0+(1:NQ),(((PARAMNUM*OUTMEASURES)-OUTMEASURES)+MEASURE),CURVENUM] +
				#		RESULTSARRAY[N0-(1:NQ),(((PARAMNUM*OUTMEASURES)-OUTMEASURES)+MEASURE),CURVENUM]
				
				#Y_VECM = RESULTSARRAY[N0+(1:NQ),(((PARAMNUM*OUTMEASURES)-OUTMEASURES)+MEASURE),CURVENUM] -
				#		RESULTSARRAY[N0-(1:NQ),(((PARAMNUM*OUTMEASURES)-OUTMEASURES)+MEASURE),CURVENUM]
			
				AC<-array(0,dim=c(1,4,1))
				BC<-array(0,dim=c(1,4,1))
				
				COMPL<-0

				rangeJ<-OMi/2
				for(j in 1:rangeJ)
				{
                    			ANGLE<-(j*2*(1:NQ)*pi/N)
                    			C_VEC<-cos(ANGLE)
                    			S_VEC<-sin(ANGLE)
					AC[j]<-(MEASURE_RESULTS_FOR_PARAM[N0] + t(Y_VECP) %*% C_VEC)/N
                    			#AC[j]<-(RESULTSARRAY[N0,(((PARAMNUM*OUTMEASURES)-OUTMEASURES)+MEASURE),CURVENUM]+t(Y_VECP) %*% C_VEC)/N
                    			BC[j] = t(Y_VECM) %*% S_VEC/N
                    			COMPL<-COMPL+AC[j]^2+BC[j]^2
                		}

				# Computation of V_{(ci)}.
               			Vci[CURVENUM]<-2*COMPL

		                # Fourier coeff. at [P*OMi, for P=1:MI].
                		COMPL = 0

				Y_VECP <- MEASURE_RESULTS_FOR_PARAM[N0+(1:NQ)] + MEASURE_RESULTS_FOR_PARAM[N0-(1:NQ)]

				Y_VECM <- MEASURE_RESULTS_FOR_PARAM[N0+(1:NQ)] - MEASURE_RESULTS_FOR_PARAM[N0-(1:NQ)]
				#Y_VECP = RESULTSARRAY[N0+(1:NQ),(((PARAMNUM*OUTMEASURES)-OUTMEASURES)+MEASURE),CURVENUM] +
						#RESULTSARRAY[N0-(1:NQ),(((PARAMNUM*OUTMEASURES)-OUTMEASURES)+MEASURE),CURVENUM]

				#Y_VECM = RESULTSARRAY[N0+(1:NQ),(((PARAMNUM*OUTMEASURES)-OUTMEASURES)+MEASURE),CURVENUM] -
						#RESULTSARRAY[N0-(1:NQ),(((PARAMNUM*OUTMEASURES)-OUTMEASURES)+MEASURE),CURVENUM]
	
				
				for(i in seq(OMi,OMi*MI,OMi))
				{
					ANGLE<-i*2*(1:NQ)*pi/N
	   	                        C_VEC<-cos(ANGLE)
                   			S_VEC<-sin(ANGLE)
					AC[j]<-(MEASURE_RESULTS_FOR_PARAM[N0] + t(Y_VECP) %*% C_VEC)/N
					#AC[j]<-(RESULTSARRAY[N0,(((PARAMNUM*OUTMEASURES)-OUTMEASURES)+MEASURE),CURVENUM]+t(Y_VECP) %*% C_VEC)/N
      		                        BC[j] = t(Y_VECM) %*% S_VEC/N
                                        COMPL<-COMPL+AC[j]^2+BC[j]^2
				}
			
				# Computation of V_i.
                		Vi[CURVENUM]<-2*COMPL
               			# AVi = AVi+Vi;
                		# Computation of the total variance in the time domain.
				V[CURVENUM]<- t(MEASURE_RESULTS_FOR_PARAM) %*% MEASURE_RESULTS_FOR_PARAM/N
				#V[CURVENUM]<- t(RESULTSARRAY[,(((PARAMNUM*OUTMEASURES)-OUTMEASURES)+MEASURE),CURVENUM]) %*% 
				#	RESULTSARRAY[,(((PARAMNUM*OUTMEASURES)-OUTMEASURES)+MEASURE),CURVENUM]/N

               		}	# END CURVE NUMBER LOOP

			# CALCULATE SENSITIVITY INDEXES
			Si[PARAMNUM,1,MEASURE]<-mean(Vi)/mean(V)
			STi[PARAMNUM,1,MEASURE]<-1-mean(Vci)/mean(V)
			rangeSi[PARAMNUM,,MEASURE]<-Vi/V 
			rangeSTi[PARAMNUM,,MEASURE]<-1-(Vci/V)      
	
			if(is.nan(Si[PARAMNUM,1,MEASURE]))
			{
				Si[PARAMNUM,1,MEASURE]<-0
			}
			if(is.nan(STi[PARAMNUM,1,MEASURE]))
			{
				STi[PARAMNUM,1,MEASURE]<-0
			}	
			for(i in seq(1:length(rangeSi[PARAMNUM,,MEASURE])))
			{
				if(is.nan(rangeSi[PARAMNUM,i,MEASURE]))
				{
					rangeSi[PARAMNUM,i,MEASURE]<-0
				}
			}	
			for(i in seq(1:length(rangeSTi[PARAMNUM,,MEASURE])))
			{
				if(is.nan(rangeSTi[PARAMNUM,i,MEASURE]))
				{
					rangeSTi[PARAMNUM,i,MEASURE]<-0
				}
			}		
			#if(is.nan(rangeSi[PARAMNUM,,MEASURE]))
			#{
			#	rangeSi[PARAMNUM,,MEASURE]<-0
			#}
			#if(is.nan(rangeSTi[PARAMNUM,,MEASURE]))
			#{
			#	rangeSTi[PARAMNUM,,MEASURE]<-0
			#}
			
        
			         

		} # END PARAMNUM
	}

	# THE Si, STi, RANGESi, AND RANGESTi ARRAYS ARE RETURNED AS A LIST
	return(list(Si=Si,STi=STi,rangeSi=rangeSi,rangeSTi=rangeSTi))
}

