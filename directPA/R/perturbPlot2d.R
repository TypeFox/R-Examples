#' Perturbation Plot
#'
#' This function takes in a matrix of test statistics with two columns (2-dimensional space) and the 
#' annotation list such as pathway annotation or kinase-substrate annotation, and visualize the enrichment
#' of pathways or kinases in direction specific manner.
#' 
#' @usage perturbPlot2d(Tc, annotation, minSize=5, ...)
#' @param Tc a numeric matrix. The columns are genes or phosphorylation sites and the columns are treatments 
#' vs control statistics.
#' @param annotation a list with names correspond to pathways or kinases and elements correspond to genes or
#' substrates belong to each pathway or kinase, respectively.
#' @param minSize the size of annotation groups to be considered for calculating enrichment. Groups 
#' that are smaller than the minSize will be removed from the analysis.
#' @param ... parameters for controling the plot.
#' @return a list of coordinates for pathways or kinases
#' @export 
#' @examples
#' # load the phosphoproteomics dataset
#' data(HEK)
#' 
#' # load the kinase-substrate annoations
#' data(PhosphoSite)
#' 
#' perturbPlot2d(Tc=HEK, annotation=PhosphoSite.mouse, cex=3)
#' 
perturbPlot2d <- function(Tc, annotation, minSize=5, ...) {

  # step 1. convert statistics into z-scores
  Tc.zscores <- apply(Tc, 2, function(x){qnorm(rank(x)/(nrow(Tc)+1))})
  
  # step 2. filter the groups that are smaller than the minimun cutoff
  DE = lapply(annotation, function(x){
    if(sum(rownames(Tc.zscores) %in% x) >= minSize) {
      X <- Tc.zscores[rownames(Tc.zscores)%in%x,]
      n = nrow(X)
		  Z1 = sum(X[,1])/sqrt(n)
		  Z2 = sum(X[,2])/sqrt(n)
		  list(Z1=Z1, Z2=Z2)
    }
  })

	# step3. filter DE that has 0 element
   DE <- DE[which(sapply(DE, length) != 0)]
   Z1 <- unlist(sapply(DE, function(x){x[1]}))
   Z2 <- unlist(sapply(DE, function(x){x[2]}))
   
   # visualization
   plot(Z1, Z2, col="darkblue", pch=16, xlab=colnames(Tc)[1], ylab=colnames(Tc)[2], ...)
   textxy(Z1,Z2, names(DE), col="black", cex=1)
   abline(v=0, h=0, col="gold", lty=2)
   abline(a=0, b=1, col="darkgreen", lty=2)
   abline(a=0, b=-1, col="darkgreen", lty=2)

   r <- ceiling(max(sqrt(Z1^2 + Z2^2)))
   for(i in seq(0, r, r/5)) {
      theta = seq(-3.14,3.14,0.05)
      lines(i*cos(theta),i*sin(theta),col = 'gray', type="l")
   }

   points(Z1, Z2, col="darkblue", pch=16, cex=2)

   ## return the results
   result <- list()
   result$Z1 <- Z1
   result$Z2 <- Z2
   return(result)
}
