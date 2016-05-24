

###########################################################
# NMI inference, helper function
BIFIE_NMI_inference_parameters <- function( parsM , parsrepM , fayfac ,
	RR , Nimp , Nimp_NMI , comp_cov = FALSE ){
	
    res0 <- .Call( "bifie_comp_vcov_within" , parsM , parsrepM , fayfac , 
				RR , Nimp , package="BIFIEsurvey" )
	u_diag <- res0$u_diag

	NV <- length(u_diag) / Nimp 
	u_diag <- array( u_diag , dim = c( NV , Nimp_NMI[2] , Nimp_NMI[1] ) )
	u_diag <- aperm( u_diag , c(3,2,1) )
	u <- array( 0 , dim=c( Nimp_NMI[1] , Nimp_NMI[2] , NV , NV ) )

	
	for (ii in seq( 1 , Nimp_NMI[1] ) ){
		for (jj in seq( 1 , Nimp_NMI[2] ) ){
			  if (NV>1){	
					diag(u[ii,jj,,]) <- u_diag[ii,jj,]		
						}
			   if (NV==1){
					u[ii,jj,1,1] <- u_diag[ii,jj,1]
							}
										}
								}								
	qhat <- array( parsM , dim=c(NV , Nimp_NMI[2] , Nimp_NMI[1] ) )
	qhat <- aperm( qhat , c(3,2,1) )								
	# inference using miceadds package
	res1 <- miceadds::NMIcombine( qhat=qhat , u = u , comp_cov= comp_cov , is_list=FALSE )	
	
	res1$df <- ifelse( res1$df > 1000 , Inf , res1$df )
	
	# output management
	res1$pars <- res1$qbar
	res1$pars_se <- sqrt( diag( res1$Tm ) )	
	res1$pars_fmi <- res1$lambda
	res1$pars_fmiB <- res1$lambda_Between
	res1$pars_fmiW <- res1$lambda_Within		
	res1$pars_varBetween1 <- diag(res1$Bm)
	res1$pars_varBetween2 <- diag(res1$Wm)
	res1$pars_varWithin <- diag(res1$ubar)	
	
	return(res1)
		}
#################################################################		