
plot.immer_HRM <- function( x , ... ){
	class(x) <- "mcmc.sirt"
	graphics::plot( x , ... )
		}
	
