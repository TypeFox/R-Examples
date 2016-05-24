#' Boxplot and ROC Curves
#'
#' Function to display boxplots and ROC curves to show model fit in terms of in-sample link prediction.
#'
#' @param Y (\code{N} x \code{N}) binary adjacency matrix, or list containing the adjacency matrices.
#' @param EZ (\code{N} x \code{D}) matrix (or list of matrices) containing the posterior means of the latent positions
#' @param xiT vector of posterior means of the parameter \eqn{\alpha}
#' @param BOXPLOT logical; if \code{TRUE} draws the boxplot. Default \code{BOXPLOT = FALSE}
#' @param ROC logical; if \code{TRUE} draws the ROC curve. Default \code{ROC = FALSE}
#' @param Lroc number of intervals in the ROC curve. Default \code{Lroc = 100}
#' @param labelsPlot main title for the boxplot. Default \code{labelsPlot = NULL}
#' @param powdist vector of power of the distance default \code{powdist = 2}, squared Euclidean distance, the alternative is 1, for the Euclidean distance
#' @param cexRocLeg \code{cex} for the ROC curve. Default \code{cexRocLeg = .8}
#' @param colRoc \code{col} for the ROC curve. Default \code{colRoc = seq(2, Ndata + 1)}
#' @param ltyRoc \code{lty} for the ROC curve. Default \code{ltyRoc = seq(2, Ndata + 1)}
#' @param lwdRoc \code{lwd} for the ROC curve. Default \code{lwdRoc = 2}
#' @param ... Arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}). 
#' @return The area under the ROC curve (AUC) and the selected plots. The closer the AUC takes values to 1 the better the fit.
#' @seealso \code{\link{lsm}}, \code{\link{lsjm}}

#' @references Gollini, I., and Murphy, T. B. (2014), "Joint Modelling of Multiple Network Views", Journal of Computational and Graphical Statistics \url{http://arxiv.org/abs/1301.3759}.
#' @export

#' @examples
#' 
#' N <- 20
#' Y <- network(N, directed = FALSE)[,]
#'
#' modLSM <- lsm(Y, D = 2) 
#' bp <- boxroc(Y, 
#'	EZ = modLSM$lsmEZ,
#'	xiT = modLSM$xiT, 
#'	Lroc = 150, 
#'	ROC = TRUE, 
#'	BOXPLOT = TRUE)
#' 
#' print(bp)

boxroc <- function(Y, EZ, xiT, BOXPLOT = FALSE, ROC = FALSE, Lroc = 100, labelsPlot = NULL, 
                   powdist = 2, cexRocLeg = .8, 	colRoc = seq(2, Ndata + 1), 
                   ltyRoc = seq(2, Ndata + 1), lwdRoc = 2, ...)
{
	if(is.matrix(Y)){ 
		Yo <- Y
		Y <- list()
		for(nd in 1:length(xiT)) Y[[nd]] <- Yo
		}
	
	stopifnot(is.list(Y), sapply(Y, is.adjacency))
	stopifnot(length(xiT) == length(Y))

	N <- nrow(Y[[1]])
	Ndata <- length(xiT)
	
	if(is.null(labelsPlot)) {
		if(Ndata > 1) labelsPlot <- paste('Network', 1:Ndata)
	}
	
	if(length(Y) == 1 & Ndata > 1){
		for(i in 2:Ndata){
			Y[[i]] <- Y[[1]]
		}	
	}	
	
	stopifnot(length(Y) == Ndata)
	
	if(length(powdist) == 1) powdist <- rep(powdist, Ndata)
	
	if(is.matrix(EZ)) EZ <- list(EZ)
	
	stopifnot(is.list(EZ))
	
	p1mod <- list()
	statbplot <- list()
	
	f0 <- exp(xiT) / ( 1 + exp(xiT) )

	for(nd in 1L:Ndata){
		num <- exp(xiT[nd] - dist(EZ[[nd]])^powdist[nd])
		p1mod[[nd]] <- as.matrix(num / (1 + num))
	}

# # # # # # BOXPLOT  # # # # # #
	if(BOXPLOT){	
			
		if(ROC){ 
			par(mfrow=c(1, Ndata + 1))
			} else {
			par(mfrow=c(1, Ndata))
			}
		
		}
		
		for(nd in 1L:Ndata){
		
		statbplot[[nd]] <- boxplot(p1mod[[nd]][Y[[nd]] == 0], p1mod[[nd]][Y[[nd]] == 1], 
					names = c(expression(italic(y[ij])==0), expression(italic(y[ij])==1)), ylim=c(0,max(f0)),
					main = labelsPlot[nd] ,
					ylab =expression(italic(p)(italic(y[ij])== paste(1, ' | Model'))), 
					plot = BOXPLOT)$stats
		
		if(BOXPLOT) abline(h = exp(xiT[nd]) / (1 + exp(xiT[nd])), col = 'red', lty = 2)
		
		}
		
		# # # # # # ROC # # # # # #
		
		trueposRate <- matrix(NA, Lroc, Ndata)
		falseposRate <- matrix(NA, Lroc, Ndata)
		
		i<-0
		
		sY <- sapply(Y, sum)
		
		Y1 <- lapply(Y, function(y) 1 - y)
	
		yp0 <- list()
		yp1 <- list()
		
		for(nd in 1L:Ndata){
			diag(Y1[[nd]]) <- 0
			yp0[[nd]] <- Y[[nd]] * p1mod[[nd]]
			yp1[[nd]] <- Y1[[nd]] * p1mod[[nd]]
			}
		
		s1YN <- sapply(Y1, sum) 
	
		for(l in seq(0, 1, length = Lroc))
		{
			i<-i+1
			
				trueposRate[i,] <- colSums(sapply(yp0, function(y) y > l)) / sY
				falseposRate[i,] <- colSums(sapply(yp1, function(y) y > l))  / s1YN
			
		}
		
		AUC <- colSums(-apply(falseposRate, 2, diff) * (apply(trueposRate, 2, diff) + 2 * trueposRate[-Lroc,]) / 2)
		
		if(ROC){ 
			
			matplot(falseposRate, trueposRate, 
			xlim = c(0,1), ylim = c(0,1), 
			type = 'l',
			main = 'ROC curve',
			xlab = 'false positive rate',
			ylab = 'true positive rate', 
			col = colRoc, lty = ltyRoc, pch = 2, lwd = 3)

			legend("bottomright",
				paste('AUC', labelsPlot, '=', round(AUC, 2) * 100,'%'),
				cex = cexRocLeg, col = colRoc, lty = ltyRoc, lwd = lwdRoc)	
			}
	
	names(statbplot) <- labelsPlot
	names(AUC) <- labelsPlot
	print(AUC)
	return(list(AUC = AUC, statsboxplot = statbplot))
}
