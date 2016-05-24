Measure.R2S <-
function(biloci,struc,na.presence=TRUE){
	CALC	=	TRUE
	M.r2s	=	NA
	
	if (na.presence==TRUE){
		if (any(is.na(biloci)==TRUE)){
			ligne	=	na.action(na.omit(biloci))
			#less than 5 non-missing data
			CALC=( length(ligne)<=(nrow(biloci)-5) )
			if(CALC==TRUE){
				biloci	=	biloci[-ligne,]
				struc	=	struc[-ligne,]
			}
		}
	}
	
	if(CALC==TRUE){
		sig_biloci			=	var(biloci)
		sig_struc			=	var(struc)
		sig_struc_inv		=   Inv.proj.matrix.sdp(sig_struc)
		sig_biloci_struc	=	cov(biloci,struc)
		sig_struc_biloci	=	t(sig_biloci_struc)
		SIG					=	sig_biloci-(sig_biloci_struc%*%sig_struc_inv%*%sig_struc_biloci)
		ifelse ((SIG[1,1]<0.0000001)|(SIG[2,2]<0.0000001),M.r2s<-0,M.r2s<-(SIG[1,2])^2/(SIG[1,1]*SIG[2,2]) )
		as.numeric(M.r2s)
	}
	
	M.r2s
}
