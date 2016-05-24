
###############################################
# summary function for mira.nmi object
summary.mira.nmi <- function(object, ...) {	
	B <- object$Nimp["between"]
	W <- object$Nimp["within"]		
	for (bb in 1:B){	
		# code adapted from summary.mira form mice package
		for (i in 1:W){
			cat("\n", "## summary of imputation: between", 
					bb , "- within" , i, ":\n")
			print(summary(object$analyses[[bb]][[i]], ...), ...)
				}
			}
}
#################################################