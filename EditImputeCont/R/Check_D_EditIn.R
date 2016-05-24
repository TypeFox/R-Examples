.Check_D_EditIn <- function(Y_temp){
	
	if ( sum(Y_temp<0) > 0 ){
		message("Negative values in the edit-in dataset are considered as missing values.")
	}
	
	if ( sum(is.na(Y_temp)) > 0 ){
		Y_temp[is.na(Y_temp)]=-999
		message("NA in the edit-in dataset are considered as missing values.")
	}
			
	if ( sum(Y_temp==0)>0 ){
		
		message("There are zeros in some variables. For compuational purpose, they are replaced with some small postive values.")
	
		for (i_var in 1:dim(Y_temp)[[2]]){
		
				whichZeros = which( Y_temp[,i_var]==0 )
				if ( length(whichZeros) > 0 ){
					Y_temp[whichZeros,i_var] = 0.01					
					if ( length(whichZeros) == 1 ){
						message( paste("Single zero in ",varnames = dimnames(Y_temp)[[2]][i_var]," is replaced by 0.01",sep="" ) )
					} else {
						message( paste(length(whichZeros)," zeros in ",varnames = dimnames(Y_temp)[[2]][i_var]," are replaced by 0.01",sep="" ) )
					}
				} # if
		
		} # for
	
	} 
		
	return(Y_temp)
	
} # Check_D_obs