#' Goodness-of-Fit diagnostics for LSM model
#'
#' This function produces goodness-of-fit diagnostics for LSM model.
#'
#' @param object object of class \code{'lsm'}
#' @param Y (\code{N} x \code{N}) binary adjacency matrix
#' @param Ysim list containing simulated (\code{N} x \code{N}) adjacency marices.  Default \code{Ysim = NULL}
#' @param nsim number of simulations. Default \code{nsim = 100}
#' @param seed for simulations 
#' @param directed if the network is directed or not.  Default \code{directed = NULL} 
#' @param stats statistics used. Default \code{stats = NULL}
#' @param doplot draw boxplot.  Default \code{doplot = TRUE}
#' @param parm do all the plots in one window. Default \code{parm = TRUE}
#' @seealso \code{\link{lsm}}, \code{\link{simulateLSM}}, \code{\link{plot.gofobj}}, \code{\link{print.gofobj}}
#' @export
#' @examples
#' 
#' Y <- network(15, directed = FALSE)[,]
#'
#' modLSM <- lsm(Y, D = 2) 
#' myGof <- goflsm(modLSM, Y = Y)

goflsm <- function(object, Y, Ysim = NULL, nsim = 100, seed, directed = NULL, stats = NULL, doplot = TRUE, parm = TRUE) {
	
	stopifnot(inherits(object, 'lsm'))
	stopifnot(is.adjacency(Y))
	stopifnot(is.logical(doplot) & length(doplot) == 1)
	stopifnot(is.logical(parm) & length(parm) == 1)
	
	nsim <- as.integer(nsim)
	stopifnot(nsim > 0)
	
	if(is.null(directed)){ 
			directed <- !all(Y == t(Y))
		} else {
			stopifnot(is.logical(directed) & length(directed) == 1)
		}
	
	if(is.null(Ysim)){
			Ysim <- simulateLSM(object = object, Y = Y, nsim = nsim, seed = seed, directed = directed)
		} else {
      is.list(Ysim)
				stopifnot(lapply(Ysim, is.adjacency))
		}
		
		nsim <- length(Ysim)
		
		if(directed){
			allstats <- c('idegree', 'odegree',  'esp', 'dsp', 'triadcensus', 'distance')	
		} else {
			allstats <- c('degree', 'esp', 'dsp', 'triadcensus', 'distance')
		}
		
		if(is.null(stats)){
			GOFstats <- allstats
		} else {
			
			if(!all(stats %in% allstats)) stop(paste("stats must be chosen between these values:", paste(allstats, collapse = ', '), sep = '\n'))

			GOFstats <- stats[stats %in% allstats]	
		}	
		
		ngofst <- length(GOFstats)
		N <- nrow(Y)
		
		SimStats <- list()
		summaryStats <- list()
	
	if( 'degree' %in% GOFstats){
		
		GOFmat <- sapply(Ysim, degree2)
		GOFmat <- cbind(GOFmat, degree2(Y))
		rownames(GOFmat) <- seq(0, N-1)
		
		SimStats[['degree']] <- GOFmat[!(rowSums(GOFmat) == 0),]
		
		obs.model <- SimStats[['degree']][, nsim + 1]
		sim.model <- SimStats[['degree']][,-(nsim + 1)]
		pval.leq <- apply(sim.model <= obs.model, 1, mean)
		pval.geq <- apply(sim.model >= obs.model, 1, mean)
		pval <- pmin(1, 2 * pmin(pval.leq, pval.geq))
		
		summaryStats[['degree']] <-  cbind(SimStats[['degree']][, nsim + 1],  
		apply(SimStats[['degree']][, -(nsim + 1)], 1, min), 
		apply(SimStats[['degree']][, -(nsim + 1)],1, mean),
		apply(SimStats[['degree']][, -(nsim + 1)], 1, max), pval)
		
		colnames(summaryStats[['degree']]) <-  c('obs', 'min', 'mean',  'max', 'MC p-value')
		
		}
		
	if( 'idegree' %in% GOFstats){
		
		GOFmat <- sapply(Ysim, indegree2)
		GOFmat <- cbind(GOFmat, indegree2(Y))
		rownames(GOFmat) <- seq(0, N-1)
		
		SimStats[['idegree']] <- GOFmat[!(rowSums(GOFmat) == 0),]
		
		obs.model <- SimStats[['idegree']][,nsim + 1]
		sim.model <- SimStats[['idegree']][,-(nsim + 1)]
		pval.leq <- apply(sim.model <= obs.model, 1, mean)
		pval.geq <- apply(sim.model >= obs.model, 1, mean)
		pval <- pmin(1, 2 * pmin(pval.leq, pval.geq))
		
		summaryStats[['idegree']] <-  cbind(SimStats[['idegree']][ ,nsim + 1],  
		apply(SimStats[['idegree']][, -(nsim + 1)], 1, min), 
		apply(SimStats[['idegree']][, -(nsim + 1)], 1, mean),
		apply(SimStats[['idegree']][, -(nsim + 1)], 1, max), pval)
		
		colnames(summaryStats[['idegree']]) <-  c('obs', 'min', 'mean',  'max', 'MC p-value')
		
		}
		
		if( 'odegree' %in% GOFstats){
		
		GOFmat <- sapply(Ysim, outdegree2)
		GOFmat <- cbind(GOFmat, outdegree2(Y))
		rownames(GOFmat) <- seq(0, N-1)
		
		SimStats[['odegree']] <- GOFmat[!(rowSums(GOFmat) == 0),]
		
		obs.model <- SimStats[['odegree']][,nsim + 1]
		sim.model <- SimStats[['odegree']][,-(nsim + 1)]
		pval.leq <- apply(sim.model <= obs.model, 1, mean)
		pval.geq <- apply(sim.model >= obs.model, 1, mean)
		pval <- pmin(1, 2 * pmin(pval.leq, pval.geq))
		
		summaryStats[['odegree']] <-  cbind(SimStats[['odegree']][, nsim + 1],  
		apply(SimStats[['odegree']][, -(nsim + 1)], 1, min), 
		apply(SimStats[['odegree']][, -(nsim + 1)], 1, mean),
		apply(SimStats[['odegree']][, -(nsim + 1)], 1, max), pval)
		
		colnames(summaryStats[['odegree']]) <-  c('obs', 'min', 'mean',  'max', 'MC p-value')
		
		}
		
		
		if( 'esp' %in% GOFstats){
		
		GOFmat <- sapply(Ysim, function(y) summary.statistics(y ~ esp(0:(N-1)), directed = directed))
		GOFmat <- cbind(GOFmat, summary.statistics(Y ~ esp(0:(N-1)), directed = directed))  / ( 2 - directed)
		rownames(GOFmat) <- seq(0, N-1)
		
		SimStats[['esp']] <- GOFmat[!(rowSums(GOFmat) == 0),]
		
		obs.model <- SimStats[['esp']][, nsim + 1]
		sim.model <- SimStats[['esp']][, -(nsim + 1)]
		pval.leq <- apply(sim.model <= obs.model, 1, mean)
		pval.geq <- apply(sim.model >= obs.model, 1, mean)
		pval <- pmin(1, 2 * pmin(pval.leq, pval.geq))
		
		summaryStats[['esp']] <-  cbind(SimStats[['esp']][, nsim + 1],  
		apply(SimStats[['esp']][, -(nsim + 1)], 1, min), 
		apply(SimStats[['esp']][, -(nsim + 1)], 1, mean),
		apply(SimStats[['esp']][, -(nsim + 1)], 1, max), pval)
		
		colnames(summaryStats[['esp']]) <-  c('obs', 'min', 'mean',  'max', 'MC p-value')
		
		}
		
		
		if( 'dsp' %in% GOFstats){
		
		GOFmat <- sapply(Ysim, function(y) summary.statistics(y ~ dsp(0:(N-1)), directed = directed))
		GOFmat <- cbind(GOFmat, summary.statistics(Y ~ dsp(0:(N-1)), directed = directed)) / (2 - directed)
		rownames(GOFmat) <- seq(0, N-1)
		
		SimStats[['dsp']] <- GOFmat[!(rowSums(GOFmat) == 0), ]
		
		obs.model <- SimStats[['dsp']][, nsim + 1]
		sim.model <- SimStats[['dsp']][, -(nsim + 1)]
		pval.leq <- apply(sim.model <= obs.model, 1, mean)
		pval.geq <- apply(sim.model >= obs.model, 1, mean)
		pval <- pmin(1, 2 * pmin(pval.leq, pval.geq))
		
		summaryStats[['dsp']] <-  cbind(SimStats[['dsp']][, nsim + 1],  
		apply(SimStats[['dsp']][-(nsim + 1),], 1, min), 
		apply(SimStats[['dsp']][-(nsim + 1),], 1, mean),
		apply(SimStats[['dsp']][-(nsim + 1),], 1, max), pval)
		
		colnames(summaryStats[['dsp']]) <- c('obs', 'min', 'mean',  'max', 'MC p-value')
		
		}
		
				
		if ("triadcensus" %in% GOFstats) {

        if(directed) {
        	
        	 GOFmat <- sapply(Ysim, function(y) summary.statistics(y ~triadcensus(0:15), directed = directed))
             GOFmat <- cbind(GOFmat, summary.statistics(Y ~ triadcensus(0:15), directed = directed))    
             rownames(GOFmat) <- c("003", "012", "102", "021D", "021U", "021C", "111D", 
             "111U", "030T", "030C",  "201", "120D", "120U", "120C", "210", "300")    
        } else {
        	
             GOFmat <- sapply(Ysim, function(y) summary.statistics(y ~triadcensus(c(0,2,10,15)), directed = directed))   
             GOFmat <- cbind(GOFmat, summary.statistics(Y ~ triadcensus(c(0,2,10,15)), directed = directed))    
             rownames(GOFmat) <- c("0", "1", "2", "3")  
        
        }

        SimStats[['triadcensus']] <- GOFmat#[,!(colSums(GOFmat) == 0)]
        
        obs.model <-SimStats[['triadcensus']][, nsim + 1]
		sim.model <- SimStats[['triadcensus']][, -(nsim + 1)]
		pval.leq <- apply(sim.model <= obs.model, 1, mean)
		pval.geq <- apply(sim.model >= obs.model, 1, mean)
		pval <- pmin(1, 2 * pmin(pval.leq,pval.geq))
		
		summaryStats[['triadcensus']] <-  cbind(SimStats[['triadcensus']][, nsim + 1],  
		apply(SimStats[['triadcensus']][, -(nsim + 1)], 1, min), 
		apply(SimStats[['triadcensus']][, -(nsim + 1)], 1, mean),
		apply(SimStats[['triadcensus']][, -(nsim + 1)], 1, max), pval)
		
		colnames(summaryStats[['triadcensus']]) <-  c('obs', 'min', 'mean',  'max', 'MC p-value')
		
		}
		
	if( 'distance' %in% GOFstats){
		
		GOFmat <- sapply(Ysim, function(y) ergm.geodistdist(as.network(y, directed = directed)))
		GOFmat <- cbind(GOFmat, ergm.geodistdist(as.network(Y, directed = directed)))
		rownames(GOFmat) <- c(seq(1, N-1), 'NR')
		
		SimStats[['distance']] <- GOFmat[!(rowSums(GOFmat) == 0), ]
		
		obs.model <- SimStats[['distance']][, nsim + 1]
		sim.model <- SimStats[['distance']][, -(nsim + 1)]
		pval.leq <- apply(sim.model <= obs.model, 1, mean)
		pval.geq <- apply(sim.model >= obs.model, 1, mean)
		pval <- pmin(1, 2 * pmin(pval.leq,pval.geq))
		
		summaryStats[['distance']] <-  cbind(SimStats[['distance']][, nsim + 1],  
		apply(SimStats[['distance']][, -(nsim + 1)], 1, min), 
		apply(SimStats[['distance']][, -(nsim + 1)], 1, mean),
		apply(SimStats[['distance']][, -(nsim + 1)], 1, max), pval)
		
		colnames(summaryStats[['distance']]) <-  c('obs', 'min', 'mean',  'max', 'MC p-value')
		
		}
		
		
		robj <- summaryStats		
		attributes(robj) <- list(SimStats = SimStats, N = N, directed = directed)
		names(robj) <- names(SimStats)
		
		list(SimStats = SimStats) 
		class(robj) <- 'gofobj'
		
	if(doplot) {
		plot(robj, parm = parm)
	}
		
	robj
}