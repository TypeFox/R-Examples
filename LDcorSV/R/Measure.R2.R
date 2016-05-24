Measure.R2 <-
function(biloci,na.presence=TRUE){
	CALC	=	TRUE
	M.r2	 =	NA
	
	if (na.presence==TRUE){
		if (any(is.na(biloci)==TRUE)){
			ligne	=	na.action(na.omit(biloci))
			#less than 5 non-missing data
			CALC=( length(ligne)<=(nrow(biloci)-5) )
			if(CALC==TRUE){
				biloci	=	biloci[-ligne,]
			}            
		}
	}
	
	if(CALC==TRUE){
		SIG=var(biloci)
		ifelse((SIG[1,1]<0.0000001) | (SIG[2,2]<0.0000001),M.r2<-0, M.r2<-(SIG[1,2])^2/(SIG[1,1]*SIG[2,2]) )
		as.numeric(M.r2)
	}
	
	M.r2
}
