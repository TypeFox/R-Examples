
MHprop_refresh <- function( MHprop ){
					vars <- MHprop$VARS_refreshing
					V <- length(vars)		
					for (vv in 1:V){
						# vv <- 1
						var.vv <- vars[vv]
						ri <- MHprop$refresh_iter[[ var.vv ]]
						accept <- MHprop$accept[[ var.vv ]] / MHprop$refresh_iter[[ var.vv ]]
						SD <- MHprop$SD[[ var.vv ]]
						SDchange <- MHprop$refresh_SDchange[[ var.vv ]]

						MHprop$SD[[ var.vv ]] <- MHprop_refresh_parstype( accept , SD , MHprop , SDchange )
						MHprop$accept[[ var.vv ]] <- 0*MHprop$accept[[var.vv]]
						MHprop$refresh_count[[var.vv]] <- 0				
									}						
				return(MHprop)
						}


##################################################################################
MHprop_refresh_parstype <- function( accept , SD , MHprop , SDchange ){
    #***********************************
	# vector
	if ( is.vector(accept) ){
		SD <- MHprop_refresh_pars( acc =accept , SD.pp=SD , MHprop ,
										SDchange )
								}
    #***********************************
	# matrix
	if ( is.matrix(accept) ){
		NP <- ncol(accept)
		for (pp in 1:NP){
			SD[,pp] <- MHprop_refresh_pars( acc =accept[,pp] , SD.pp=SD[,pp] , MHprop ,
										SDchange )
							}
						}
    #***********************************						
	return(SD)	
				}			
##################################################################################

##################################################################################							
MHprop_refresh_pars <- function( acc , SD.pp , MHprop , SDchange ){
			target <- mean( MHprop$accept_bounds )
       if (MHprop$refresh_formula){
			SD.pp <- ifelse( acc < MHprop$accept_bounds[1] ,
						SD.pp / ( 2 - acc / target ) , SD.pp )	
			SD.pp <- ifelse( acc > MHprop$accept_bounds[2] ,
						SD.pp * ( 2 - (1-acc)/(1-target) ) , SD.pp )	
						} else {						
			SD.pp <- ifelse( acc < MHprop$accept_bounds[1] , SD.pp - SDchange  , SD.pp )
			SD.pp <- ifelse( acc > MHprop$accept_bounds[2] , SD.pp + SDchange  , SD.pp )
			SD.pp <- ifelse( SD.pp < SDchange , SDchange , SD.pp )
							}
			return(SD.pp)
					}
##################################################################################					