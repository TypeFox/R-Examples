#' Create a heatmap with censored values
#'
#' Instead of using the whole range of expression values
#' we here limit the heatmap to \eqn{100 - 2* cutoff} percent of the
#' data. All values below the cutoff percentile and above 1 -
#' the cutoff percentile will be set to the same value.
#'
#' @param data a normalized but not scaled matrix
#' @param distance the distance measure to be used. One of c("eucsq", "euclidean","maximum","manhattan","canberra","binary","minkowski","correlation")
#' @param dendros which dedrograms should be computed. One of c("row","column","both","none")
#' @param cutoff the percentile(s) where to censor the data. If one value, the cutoff is used as lower percentile and the upper cutoff is set to 1-cutoff.
#'        Otherwise, if two values are provided the first is used for the lower percentile and the second for the upper percentile.
#' @param toFile should the heatmap plotted into a PDF-file
#' @param fileName a file name
#' @param cols a vector of colors
#' @param hclust.method A method to be used for hierarchical clustering. See help for hclust.
#' @param ... additional parameters passed to heatmap.2
#' @return the a list containing the row and column dendrogram of the matrix
#' @note \code{plotHeatmap} is a wrapper for heatmap.2
#' it does a row wise scaling.
#' @author Marc Johannes \email{M.Johannes@@dkfz.de}
plotHeatmap <- function(data, distance = "eucsq", dendros="both", cutoff=0.005, toFile=FALSE, fileName="Heatmap.pdf", cols=colorpanel(100, low="blue",mid="yellow",high="red"), hclust.method="ward", scale="row", ...){


  distance.measures   <- c("eucsq", "euclidean","maximum","manhattan","canberra","binary","minkowski")
  dendro.orientation  <- c("row","column","both","none")

  if(!(distance %in% c(distance.measures, "correlation"))){stop(paste("Distance ", distance, " not implemented, yet\n", sep=""))}

  if(!(dendros %in% dendro.orientation)) { stop(paste("dendros must be one of: ",paste(dendro.orientation, collapse=","), sep=""))}

  if(scale=="row") {
	## scale the rows of the data  
	data.scaled = t(scale(t(data)))
  } else if(scale=="col") {
	data.scaled = scale(data)
  } else if(scale=="both") {
	data.scaled = t(scale(t(scale(data))))
  } else {
	data.scaled = data
  }

  dd <- col.dendrogram <- row.dendrogram <- NULL

  if(dendros %in% c("column", "both")){
    ## distance by default computes the distance between
    ## the rows of data. Thus it must be transposed!
    ## For correlation distance, cor computes correlations between the 
    ## columns of a matrix, thus transposing can be omitted
    if(distance %in% distance.measures)   {
		if(distance == "eucsq") {
			dd <- dist(t(data.scaled), method = "euclidean")^2
		} else {
			dd <- dist(t(data.scaled), method = distance)
		}
	} else if(distance == "correlation") {
		dd <- as.dist(1 - cor(data.scaled))
	}
	
    col.dendrogram <- as.dendrogram(hclust(dd,method=hclust.method))
  }

  if(dendros %in% c("row", "both")){
    ## dist gets row distances, no transpose
    ## cor uses column correlations, thus transpose the data matrix
    if(distance %in% distance.measures)   {
    	if(distance == "eucsq") {
			dd <- dist(data.scaled, method = "euclidean")^2
		} else {
			dd <- dist(data.scaled, method = distance)
		}
	} else if(distance == "correlation") {
		dd <- as.dist(1 - cor(t(data.scaled)))
	}
	
    row.dendrogram <- as.dendrogram(hclust(dd, method=hclust.method))
  }

  ## censor the data, i.e. use only
  ## 99% of it and set the rest to the
  ## 99.5 and 0.05 percentile, respectively.
  lower.cutoff <- upper.cutoff <- 1
  if(length(cutoff == 1)){
    lower.cutoff <- cutoff
    upper.cutoff <- 1-cutoff
  }
  else if(length(cutoff == 2)){
    lower.cutoff <- cutoff[1]
    upper.cutoff <- cutoff[2]
  }
  else
    stop("length of cutoff must be either 1 or 2!")
                                           
  data.censored <- data.scaled
  lower.quant   <- quantile(data.censored, lower.cutoff)
  upper.quant   <- quantile(data.censored, upper.cutoff)

  print(paste("Discarding ",sum(data.censored < lower.quant), "values below lower ", lower.cutoff*100, "% quantile and ",sum(data.censored > upper.quant), " values above upper ", upper.cutoff*100, "% quantile",sep=""))
  
  data.censored[data.censored < lower.quant] <- lower.quant
  data.censored[data.censored > upper.quant] <- upper.quant

  if(toFile) pdf(file=fileName)

  heatmap.2(data.censored, Colv=col.dendrogram, Rowv=row.dendrogram, col=cols ,scale="none", symkey=FALSE, trace="none", cexRow=0.3, cexCol=0.3, dendrogram=dendros,...)
  
  if(toFile){
    dev.off()
    print(paste("Created file ", getwd(),"/",fileName,sep=""))
  }
  
  list(row.dendrogram=row.dendrogram, column.dendrogram=col.dendrogram)
}
