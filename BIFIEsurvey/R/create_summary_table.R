
##################################################################
create_summary_table <- function( res_pars , parsM , parsrepM , dfr , BIFIEobj){	

	dfr$est <- res_pars$pars
	dfr$SE <- res_pars$pars_se
	dfr$t <- round( dfr$est / dfr$SE , 2 )
	dfr$df <- rubin_calc_df( res_pars , BIFIEobj$Nimp )	
	dfr$p <- stats::pt( - abs( dfr$t ) , df=dfr$df) * 2
	# dfr$p <- pnorm( - abs( dfr$t ) ) * 2
	dfr$fmi <- res_pars$pars_fmi
	dfr$VarMI <- res_pars$pars_varBetween
	dfr$VarRep <- res_pars$pars_varWithin

	#*****
	# inference nested multiple imputation
	if ( BIFIEobj$NMI ){
		res1 <- BIFIE_NMI_inference_parameters( parsM= parsM, parsrepM= parsrepM , 
					fayfac=BIFIEobj$fayfac , RR=BIFIEobj$RR , Nimp=BIFIEobj$Nimp , 
					Nimp_NMI=BIFIEobj$Nimp_NMI , comp_cov = FALSE )	
		dfr$fmi <- dfr$cov_VarMI <- NULL	
		dfr$est <- res1$pars
		dfr$SE <- res1$pars_se
		dfr$t <- round( dfr$est / dfr$SE , 2 )
		dfr$df <- res1$df
		dfr$p <- stats::pt( - abs( dfr$t ) , df=dfr$df) * 2			
		dfr$fmi <- res1$pars_fmi
		dfr$fmi_St1 <- res1$pars_fmiB
		dfr$fmi_St2 <- res1$pars_fmiW
		dfr$VarMI_St1 <- res1$pars_varBetween1
		dfr$VarMI_St2 <- res1$pars_varBetween2
		dfr$VarRep <- res1$pars_varWithin
	    dfr$VarMI <- res1$pars_varBetween1 + res1$pars_varBetween2		
						}	
	
	
    return(dfr)
		}
##################################################################		