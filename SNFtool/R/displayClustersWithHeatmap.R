## Arguments:
## W: Similarity matrix
## group: a numeric vector containing the groups information for each sample in W such as the result of the spectralClustering function. The order should correspond to the sample order in W.
## ColSideColors:  (optional) character vector of length ncol(x) containing the color names for a horizontal side bar that may be used to annotate the columns of x, used by the heatmap function,
## OR a character matrix with number of rows matching number of rows in x.  Each column is plotted as a row similar to heatmap()'s ColSideColors by the heatmap.plus function.
## ... other paramater that can be pass on to the heatmap (if ColSideColor is a NULL or a vector) or  heatmap.plus function (if ColSideColors is matrix)

## Details:
## Using the heatmap or heatmap.plus function to display the similarity matrix
## For representation purpose, the similarity matrix diagonal is set to the median value of W, the matrix is normalised and W = W + t(W) is applied
## In this presentation no clustering method is ran the samples are ordered in function of their group label present in the group arguments.

## Values:
## Plots the similarity matrix using the heatmap function. Samples are ordered by the clusters provided by the argument groups with sample information displayed with a color bar if the ColSideColors argument is informed.
## Autors:

displayClustersWithHeatmap <- function (W, group,ColSideColors=NULL,...) {
  normalize <- function(X) X/rowSums(X)
  ind = sort(as.vector(group), index.return = TRUE)
  ind = ind$ix
  ## diag(W) = 0
  diag(W) = median(as.vector(W))
  W = normalize(W)
  W = W + t(W)
  if(is.null(ColSideColors)){
    heatmap(W[ind, ind],scale="none",Rowv=NA,Colv=NA,...)
  }
  else{
    if(is.vector(ColSideColors)){
      heatmap(W[ind, ind],scale="none",Rowv=NA,Colv=NA,ColSideColors=ColSideColors[ind],...)
    }
    else{
      heatmap.plus(W[ind, ind],scale="none",Rowv=NA,Colv=NA,ColSideColors=ColSideColors[ind,],...)
    }
  }
}