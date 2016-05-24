
##############################################################
# summary method	
summary.gdina.wald <- function(object, digits=3, 
		vars = c("X2" , "p" , "sig" , "RMSEA" , "wgtdist") , ...){
		stats <- object$stats
		cn <- colnames(stats)
		cn <- cn[-1]
		sels <- c("NAttr")
		for ( rule in object$cdm_rules){
			sels <- c( sels , paste0( rule , "_" , vars ) )
					}
		stats <- stats[ , sels ]
		cn <- colnames(stats)	
		cn <- cn[ - c( 1 , grep("_sig" , cn) ) ]
		for (vv in cn){ 
				stats[,vv] <- round(stats[,vv],digits) 
						}	
		rownames(stats) <- paste(object$stats$item)
		print(stats)
			}
###############################################################	