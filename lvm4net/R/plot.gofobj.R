#' Plot GoF object
#'
#' Function to plot an object of class \code{'gofobj'}
#'
#' @param x object of class \code{"gofobj"}
#' @param parm do all in one plots
#' @param ... other
#' @export
#' @examples
#' Y <- network(20, directed = FALSE)[,]
#'
#' modLSM <- lsm(Y, D = 2) 
#' myGof <- goflsm(modLSM, Y = Y, doplot = FALSE)
#' plot(myGof)

plot.gofobj <- function(x, parm = TRUE, ...)
{		
	
	stopifnot(inherits(x, 'gofobj'))
	stopifnot(is.logical(parm) & length(parm) == 1)
	
	SimStats <- attributes(x)$SimStats
	N <- attributes(x)$N	
	directed <- attributes(x)$directed	
	nd <- N * (N - 1) / (2 - directed)

	stats <- names(SimStats)
	nsim <- nrow(SimStats[[1]]) - 1
	ngofst <- length(stats)
	
	if(directed){
			allstats <- c('idegree', 'odegree',  'esp', 'dsp', 'triadcensus', 'distance')	
			allxlab <- c('in degree', 'out degree',  'edge-wise shared partners', 'dyad-wise shared partners', 'triad census', 'minimum geodesic distance')	
			allylab <- c('proportion of nodes', 'proportion of nodes', 'proportion of edges', 'proportion of dyads', 'proportion of triads', 'proportion of dyads')
			
	} else {
			allstats <- c('degree', 'esp', 'dsp', 'triadcensus', 'distance')
			allxlab <- c('degree',  'edge-wise shared partners', 'dyad-wise shared partners', 'triad census', 'minimum geodesic distance')	
			allylab <- c('proportion of nodes', 'proportion of edges', 'proportion of dyads', 'proportion of triads', 'proportion of dyads')
	}			
			names(allxlab) <- allstats
			names(allylab) <- allstats

			GOFxlab <- allxlab[stats %in% allstats]	
			GOFylab <- allylab[stats %in% allstats]		

	sumstats <- list()
	
	if(parm){
		if(ngofst == 2 ){
			par(mfrow = c(1,2))
		}
		
		if(ngofst %in% c(3,4) ){
			par(mfrow = c(2,2))
		}
		
		if(ngofst >  4 ){
			par(mfrow = c(2,3))
		}
	}
	
	for(gf in stats) {
		
		if(gf != 'distance'){
	
			if(gf %in% c('degree', 'idegree', 'odegree')){
		
			sumstats <-  cbind(x[[gf]][,1],  
				apply(SimStats[[gf]][, -(nsim + 1)], 1, quantile, probs = .025), 
				apply(SimStats[[gf]][, -(nsim + 1)], 1, mean),
				apply(SimStats[[gf]][, -(nsim + 1)], 1, quantile, probs = .975)) / N
	
			simstats <- t(SimStats[[gf]][, -(nsim + 1)]) / N
			}
		
		if(gf == 'esp'){
		
		ne <- colSums(attributes(x)$SimStats[['esp']])
		
		sumstats <- cbind(x[[gf]][,1] / ne[nsim + 1],  
			apply(t(SimStats[[gf]][, -(nsim + 1)]) / ne[ -(nsim + 1)], 2, quantile, probs = .025), 
			apply(t(SimStats[[gf]][, -(nsim + 1)]) / ne[ -(nsim + 1)], 2, mean),
			apply(t(SimStats[[gf]][, -(nsim + 1)]) / ne[ -(nsim + 1)], 2, quantile, probs = .975))		
	
		simstats <- t(SimStats[[gf]][, -(nsim + 1)]) / ne[ -(nsim + 1)]
		}
		
	if(gf == 'dsp'){
			
		sumstats <-  cbind(x[[gf]][,1],  
			apply(SimStats[[gf]][, -(nsim + 1)], 1, quantile, probs = .025), 
			apply(SimStats[[gf]][, -(nsim + 1)], 1, mean),
			apply(SimStats[[gf]][, -(nsim + 1)], 1, quantile, probs = .975)) / nd
	
		simstats <- t(SimStats[[gf]][, -(nsim + 1)]) / nd
		}
		
		
	if(gf == 'triadcensus'){	
		
		nt <- sum(x[['triadcensus']][,1])
		sumstats <-  cbind(x[[gf]][,1],  
			apply(SimStats[[gf]][, -(nsim + 1)], 1, quantile, probs = .025), 
			apply(SimStats[[gf]][, -(nsim + 1)], 1, mean),
			apply(SimStats[[gf]][, -(nsim + 1)], 1, quantile, probs = .975)) / nt
	
		simstats <- t(SimStats[[gf]][, -(nsim + 1)]) / nt
		
		}		
		
		boxplot(simstats, xlab = GOFxlab[gf] , ylab = GOFylab[gf],  
				ylim = range(simstats, sumstats[,1]), outline = FALSE) 
		
		matlines(sumstats[, c(2, 4)], col = 'gray', lty = 1)
		matpoints(sumstats[, c(2,4)], col = 1, pch = 1)

		lines(sumstats[,1], col = 'red', lwd = 2)
		
		} else { ## gf == ' distance'
			
		obsstats <- x[[gf]][,1] / nd	
		simstatsCI <-  t(apply(SimStats[[gf]][, -(nsim + 1)], 1, quantile, probs = c(.025, .975))) / nd
		simstats <- t(SimStats[[gf]][, -(nsim + 1)]) / nd
			
		nssd <- length(obsstats)
			
		boxplot(simstats, xlab = GOFxlab[gf] , ylab =GOFylab[gf],  ylim = range(simstats, obsstats), outline = FALSE) 
		
		matlines(simstatsCI[-nssd, ], col = 'gray', lty = 1)
		matpoints(simstatsCI, col = 1, pch = 1)
		
		points(rep(nssd,2), simstatsCI[nssd, ] / nd, col = 'gray', pch = '_')

		lines(obsstats[-nssd], col = 'red', lwd = 2)
		points(nssd, obsstats[nssd], col = 'red', pch = 16, cex = 1)
	
		}
		
		}
		
		mtext("Goodness-of-fit diagnostics", side = 3, outer = TRUE, cex = 1.5, padj = 2)
		
		if(parm) par(mfrow = c(1,1))
		
}