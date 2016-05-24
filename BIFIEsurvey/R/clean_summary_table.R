

clean_summary_table <- function( dfr , RR , se , Nimp ){

	#****
	# RR == 0
	if ( ( ! se ) &  ( RR==0 )  ){				
	    vars <- c("df" , "cor_SE" , "cor_fmi" , "cor_VarMI" ,
					"cor_VarRep" , "cor_VarMI_St2" , "cor_VarMI_St1" ,
					"cor_fmi_St1" , "cor_fmi_St2" , "t" , "p" ,
					"cov_SE" , "cov_fmi" , "cov_VarMI" , "cov_VarRep",
					"cov_VarMI_St2" , "cov_VarMI_St1"	,
					"cov_fmi_St1" , "cov_fmi_St2" , "SE" ,
					"Var_MI" , "VarRep"
					)
		for (vv in vars){  dfr[,vv]  <- NULL }
				}
    #***
	# Nimp == 1
	if ( Nimp == 1  ){				
	    vars <- c("cov_fmi" , "cov_VarMI" , "fmi" , "VarMI")
		for (vv in vars){  dfr[,vv]  <- NULL }
				}

				
	return(dfr)
		}