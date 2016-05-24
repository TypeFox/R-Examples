efast_cvmethod <-
function(Si, rangeSi, STi, rangeSTi, OUTMEASURE,NUMPARAMS,NUMCURVES,NUMOUTMEASURES)
{
	meanSi<-array(0,dim=c(NUMPARAMS,NUMOUTMEASURES,1))
	meanSTi<-array(0,dim=c(NUMPARAMS,NUMOUTMEASURES,1))
	stdSi<-array(0,dim=c(NUMPARAMS,NUMOUTMEASURES,1))
	stdSTi<-array(0,dim=c(NUMPARAMS,NUMOUTMEASURES,1))
	errorSi<-array(0,dim=c(NUMPARAMS,1))
	errorSTi<-array(0,dim=c(NUMPARAMS,1))


	if(NUMOUTMEASURES==1)
	{
		for(PARAM in 1:NUMPARAMS)
		{
			meanSi[PARAM,OUTMEASURE,1] <- mean(rangeSi[PARAM,,OUTMEASURE])
			meanSTi[PARAM,OUTMEASURE,1] <- mean(rangeSTi[PARAM,,OUTMEASURE])
	        	
	        	stdSi[PARAM,OUTMEASURE,1] <- (sd(rangeSi[PARAM,,OUTMEASURE]))
	        	stdSTi[PARAM,OUTMEASURE,1] <- (sd(rangeSTi[PARAM,,OUTMEASURE]))

			errorSi[PARAM,1] <- sd(rangeSi[PARAM,,OUTMEASURE])/sqrt(length(rangeSi[PARAM,,OUTMEASURE]))
			errorSTi[PARAM,1] <- sd(rangeSTi[PARAM,,OUTMEASURE])/sqrt(length(rangeSTi[PARAM,,OUTMEASURE]))
	    	}
		a <- Si[,OUTMEASURE,1]/meanSi[,OUTMEASURE,1]
		b <- STi[,OUTMEASURE,1]/meanSTi[,OUTMEASURE,1]	    	
	} else
	{
		for(PARAM in 1:NUMPARAMS)
        	{
			meanSi[PARAM,OUTMEASURE,1] <- drop(mean(rangeSi[PARAM,,OUTMEASURE]))
			meanSTi[PARAM,OUTMEASURE,1] <- drop(mean(rangeSTi[PARAM,,OUTMEASURE]))
			
      			qw<-rangeSi[PARAM,,OUTMEASURE]
	   		stdy=sd(rangeSi[PARAM,,OUTMEASURE])
            		stdSi[PARAM,OUTMEASURE,1]=drop(sd(rangeSi[PARAM,,OUTMEASURE]))
            		stdSTi[PARAM,OUTMEASURE,1]=drop(sd(rangeSTi[PARAM,,OUTMEASURE]))

			errorSi[PARAM,1] <- sd(rangeSi[PARAM,,OUTMEASURE])/sqrt(length(rangeSi[PARAM,,OUTMEASURE]))
			errorSTi[PARAM,1] <- sd(rangeSTi[PARAM,,OUTMEASURE])/sqrt(length(rangeSTi[PARAM,,OUTMEASURE]))
   		}

		a <- stdSi[,OUTMEASURE,1]/meanSi[,OUTMEASURE,1]
		b <- stdSTi[,OUTMEASURE,1]/meanSTi[,OUTMEASURE,1]
	}

	CVSi <- 100*t(a)
	CVSTi <- 100*t(b)

	return(list(errorSi=errorSi,errorSTi=errorSTi,CVSi=CVSi,CVSTi=CVSTi))
}

