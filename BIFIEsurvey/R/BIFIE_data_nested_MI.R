
###########################################################
# subfunction for handling nested multiple imputation
BIFIE_data_nested_MI <- function( data.list , NMI ){
	Nimp_B <- length(data.list)
	Nimp_W <- length(data.list[[1]])
	Ntot <- Nimp_B * Nimp_W
	data.list0 <- data.list
	Nimp_NMI <- NULL
	if (NMI){
		data.list <- as.list( 1:Ntot)
		hh <- 1
		for (ii in 1:Nimp_B){
		   for (jj in 1:Nimp_W){
				data.list[[hh]] <- data.list0[[ii]][[jj]]
				hh <- hh + 1
							}
						}		
		Nimp_NMI <- c( Nimp_B , Nimp_W )
			}
	res <- list( data.list = data.list , Nimp_NMI = Nimp_NMI )
	return(res)
			}
###############################################################			
