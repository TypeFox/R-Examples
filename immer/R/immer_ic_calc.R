
immer_IC_calc <- function(ic){
    	# AIC
        ic$AIC <- ic$dev + 2*ic$np
		# AIC3
		ic$AIC3 <- ic$dev + 3*ic$np		
        # BIC
        ic$BIC <- ic$dev + ( log(ic$n) )*ic$np
		# adjusted BIC 
		ic$aBIC <- ic$dev + ( log( ( ic$n -2 ) / 24 ) )*ic$np
        # CAIC (consistent AIC)
        ic$CAIC <- ic$dev + ( log(ic$n) + 1 )*ic$np
		# corrected AIC
        ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )	
		return(ic)
				}