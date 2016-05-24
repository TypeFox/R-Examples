#' Collapse a large data set for \code{heatmap.fit}
#'
#' Reduces the size of large binary data sets by binning them according to their predicted probability [0, 1].
#'
#' @param y A vector of observations of the dependent variable (in \{0,1\}).
#' @param pred A vector of predicted Pr(y = 1) corresponding to each element of \code{y}.
#' @param init.grid The number of bins on the interval [0, 1] to use for compression of \code{pred}.
#'
#' @return A list with the elements:
#' \item{y.out}{The value of \code{y}, 0 or 1.} 
#' \item{pred.out}{The (binned) predicted Pr(y = 1) matching each observation.}
#' \item{weight.out}{A weight parameter indicating the proportion of observations in the bin; sums to one.}
#' \item{pred.total.out}{A vector of unique Pr(y = 1) bin values.}
#' \item{n.out}{The number of observations (non-empty bins) after the data are collapsed.}
#' \item{retained.obs}{A vector of indices for non-empty candidate bins (for internal use by \code{heatmap.fit}).}
#' 
#' @author Justin Esarey <justin@@justinesarey.com>
#' @export


heatmap.compress <- function(y, pred, init.grid){
   
  # cut prediction space into a smaller number of cells and classify
  # observations into these cells
  pred.init <- (cut(pred, breaks=init.grid))
  pred.idx <- as.numeric(pred.init)
  
  # extract the central location of these cells
  labs <- levels(pred.init)
  break.mat <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
        upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
  pred.out.t<-pmin(pmax(apply(X=break.mat, FUN=mean, MARGIN=1),0),1)
  
  # calculate avg. y in each cell
  tab <- table(pred.idx, y)
  y.out <- c(rep(0, length(tab[,1])), rep(1, length(tab[,2])) )
  
  pred.out <- rep( pred.out.t[as.numeric(rownames(tab))], 2)
  
  # calculate % of obs in each cell
  weight.out <- prop.table(tab)
 
  rem <- which( weight.out !=0 )
  y.out <- y.out[rem]
  pred.out <- pred.out[rem]
  weight.out <- weight.out[rem]
  
  # return compressed data set and weights
  list.out <- list(y.out = y.out, pred.out = pred.out, weight.out = weight.out, pred.total.out = sort(pred.out.t[unique(pred.idx)]), n.out = as.vector(margin.table(tab, 1)), retained.obs = rem)
  return(list.out)
}
