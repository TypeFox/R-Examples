efast_setfreq <-
function(NUMPARAMS,OMciMAX,PARAMNUM)
{
	OMci<-array(0,dim=c(1,NUMPARAMS,1))

	if(NUMPARAMS==1)
	{
    		OMci<-1
	} else if(OMciMAX==1) {
		OMci <- sample(1:1,NUMPARAMS,replace=TRUE)
	} else {
		if(OMciMAX<NUMPARAMS)
		{
			INFD <- OMciMAX			
		} else {
			INFD <- NUMPARAMS
		}
		
		ISTEP = round((OMciMAX-1)/(INFD-1))

		if(OMciMAX==1)
		{
        		ISTEP <- 0
		}

		OTMP <- 1:ISTEP:INFD*ISTEP

		fl_INFD <- floor(INFD)

		for(i in 1:NUMPARAMS)
		{
        		j = (i-1 %% fl_INFD)+1
        		OMci[i] = OTMP[j];
		}
	}
	return(OMci)
}

