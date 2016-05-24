#' Print GoF object
#'
#' Function to print an object of class \code{'gofobj'}
#'
#' @param x object of class \code{'gofobj'}
#' @param ... other
#' @export
#' @examples
#' Y <- network(20, directed = FALSE)[,]
#'
#' modLSM <- lsm(Y, D = 2) 
#' myGof <- goflsm(modLSM, Y = Y, doplot = FALSE)
#' print(myGof)

print.gofobj <- function(x,...)
{		
	
	stopifnot(inherits(x, 'gofobj'))
	
	GOFstats <- names(x)
	
	
		if('degree' %in% GOFstats){
			 cat(paste("\nGoodness-of-fit for", 'degree', '\n\n'))
			 print(round(x[['degree']], 2))
		}
		
		if('idegree' %in% GOFstats){
			cat(paste("\nGoodness-of-fit for", 'in degree', '\n\n'))
			print(round(x[['idegree']], 2))
		}
		
		if('odegree' %in% GOFstats){
			 cat(paste("\nGoodness-of-fit for", 'out degree', '\n\n'))
			 print(round(x[['odegree']], 2))
		}
		
		if('esp' %in% GOFstats){
			cat(paste("\nGoodness-of-fit for", 'edge-wise shared partners', '\n\n'))
			 print(round(x[['esp']], 2))
		}
		
		if('dsp' %in% GOFstats){
			cat(paste("\nGoodness-of-fit for", 'dyad-wise shared partners', '\n\n'))
			 print(round(x[['dsp']], 2))
		}
		
		if('triadcensus' %in% GOFstats){
			cat(paste("\nGoodness-of-fit for", 'triad census', '\n\n'))
			 print(round(x[['triadcensus']], 2))
		}
		
		if('distance' %in% GOFstats){
			cat(paste("\nGoodness-of-fit for", 'minimum geodesic distance', '\n\n'))
			print(round(x[['distance']], 2))
		}
		
}