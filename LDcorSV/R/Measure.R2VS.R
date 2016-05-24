Measure.R2VS <-
function(biloci,V,struc,na.presence=TRUE,V_inv=NULL){
	CALC	=	TRUE
	M.r2vs	=	NA
	
	if (na.presence==TRUE){
		if (any(is.na(biloci)==TRUE)){
			ligne	=	na.action(na.omit(biloci))
			#less than 5 non-missing data
			CALC=( length(ligne)<=(nrow(biloci)-5) )
			if(CALC==TRUE){
				biloci	=	biloci[-ligne,]
				struc	=	struc[-ligne,]
				V		=	V[-ligne,-ligne]
			}
		}
	}
	
	if(CALC==TRUE){
    if(is.null(V_inv)){
      V_inv  		=	Inv.proj.matrix.sdp(V) 
      rownames(V_inv)	=	rownames(V)
      colnames(V_inv)	=	colnames(V)
    }
		
		DATA			=	as.matrix(cbind(biloci,struc))
		FACT			=	(rep(1,nrow(DATA))/rep(1,nrow(DATA))%*%V_inv%*%rep(1,nrow(DATA)))%*%t(rep(1,nrow(DATA)))%*%V_inv
		MAT				=	DATA-FACT%*%DATA
		SIG				=	t(MAT)%*%V_inv%*%MAT
		
		inv_sig_struc	=	Inv.proj.matrix.sdp(SIG[-c(1,2),-c(1,2)])
		sig_l1_struc	=	SIG[1,-c(1,2)]
		sig_l2_struc	=	SIG[2,-c(1,2)]
		num_r2vs		=	(SIG[1,2]-(sig_l1_struc%*%inv_sig_struc%*%sig_l2_struc))^2
		denom11_r2vs	=	SIG[1,1]-(sig_l1_struc%*%inv_sig_struc%*%sig_l1_struc)
		denom22_r2vs	=	SIG[2,2]-(sig_l2_struc%*%inv_sig_struc%*%sig_l2_struc)
		
		ifelse ((denom11_r2vs<0.0000001)|(denom22_r2vs<0.0000001),M.r2vs<-0,M.r2vs<-num_r2vs/(denom11_r2vs*denom22_r2vs) )
		as.numeric(M.r2vs)
	}
	
	M.r2vs
}
