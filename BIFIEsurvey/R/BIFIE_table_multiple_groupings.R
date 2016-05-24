
#####################################################################
BIFIE_table_multiple_groupings <- function( dfr , res00 ){
	GR <- res00$GR
	if (GR>1){
		ind1 <- which( colnames(dfr) == "groupvar" )
		ind2 <- which( colnames(dfr) == "groupval" )
		N2 <- ncol(dfr)			
		dfr1 <- dfr[ , seq( 1 , ind1 - 1 ) , drop=FALSE ]
		for (gg in 1:GR){
			# gg <- 1
			dfr1[ , paste0("groupvar" , gg ) ] <- paste(res00$group_orig[gg])
			ind <- match( dfr$groupval , res00$group_values)
			dfr1[ , paste0("groupval" , gg ) ] <- res00$group_values_recode[ ind ,gg]		
						}
		dfr <- cbind( dfr1 , dfr[ , seq( ind2 + 1 , N2 ) , drop=FALSE] )
				}
	return(dfr)
			}
#######################################################################