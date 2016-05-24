Measure.R2V <-
function(biloci,V,na.presence=TRUE,V_inv=NULL){
	CALC	=	TRUE
	M.r2v	=	NA
	
	if (na.presence==TRUE){
		if (any(is.na(biloci)==TRUE)){
			ligne	=	na.action(na.omit(biloci))
			#less than 5 non-missing data
			CALC=( length(ligne)<=(nrow(biloci)-5) )
			if(CALC==TRUE){
				biloci	=	biloci[-ligne,]
				V		=	V[-ligne,-ligne] 
			}
		}
	}
	
	if(CALC==TRUE){
    if(is.null(V_inv)){
	  	V_inv			=	Inv.proj.matrix.sdp(V) 
		  rownames(V_inv)	=	rownames(V)
		  colnames(V_inv)	=	colnames(V)
    }
		DATA 			=	as.matrix(biloci)
		FACT			=	(rep(1,nrow(DATA))/rep(1,nrow(DATA))%*%V_inv%*%rep(1,nrow(DATA)))%*%t(rep(1,nrow(DATA)))%*%V_inv
		MAT				=	DATA-FACT%*%DATA
		SIG				=	t(MAT)%*%V_inv%*%MAT
		ifelse ((SIG[1,1]<0.0000001)|(SIG[2,2]<0.0000001),M.r2v<-0,M.r2v<-(SIG[1,2])^2/(SIG[1,1]*SIG[2,2]))
		as.numeric(M.r2v)  
	}
	
	M.r2v
}
