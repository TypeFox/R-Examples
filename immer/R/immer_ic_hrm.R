
#######################################################
# information criteria
immer_ic_hrm <- function( ic , summary.mcmcobj ){
		ic$n <- ic$N
		pars <- paste(summary.mcmcobj$parameter)		
		vars <- c("mu" , "sigma" , "a" , "b" , 
							"phi" , "psi")
		VV <- length(vars)
		Npars <- NULL
		for (vv in 1:VV){
			# vv <- 1		
			ind <- which( substring( pars , 1 , nchar( vars[vv] ) ) == vars[vv] )
			Npars[ vars[vv] ] <- length(ind)
						}
		ic$Npars <- Npars
		ic$np <- sum(Npars)
		ic <- immer_IC_calc(ic)
		return(ic)
				}
############################################################				