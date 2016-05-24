#' generate heat map for a matrix
#' 
#' \code{plotConfmat} generate heat map for a matrix or a win-loss probability matrix
#' 
#' @param conf.mat an N-by-N matrix. Either a conflict matrix or 
#' a win-loss probability matrix (the second element from \code{conductance} output)
#' @param ordering a reordering of the rows/columns, specified by a permutation of 1:N
#' @param labels if TRUE, displaying the agent names as 
#' specified in the \code{rownames()} of \code{conf.mat()} on the heatmap
#' @param ... Further argument may be supplied and processed by \code{lattice::levelplot}.
#' @return A heatmap
#' @seealso \code{\link{as.conflictmat}}, \code{\link{conductance}}
#' @examples
#' # convert an edgelist to conflict matrix
#' confmatrix <- as.conflictmat(sampleEdgelist)
#' # find win-loss probability matrix
#' perm2 <- conductance(confmatrix, 2)
#' # plotting
#' plotConfmat(perm2$p.hat)
#' @export

plotConfmat = function(conf.mat, ordering = NA, labels = FALSE, ...){
  
  # making sure input is correct
  if (!(is.matrix(conf.mat))) {
    stop("conf.mat should be a matrix.")
  }
  
  if(length(rownames(conf.mat)) == 0){
    labels = FALSE
  }
  
  if(length(ordering) == 1){
    ordering = 1:ncol(conf.mat)
  }
  
  conf.mat.ord = conf.mat[ordering, ordering]
  ramp = colorRamp(c("white","blue", "orange", "red"))
  colors = rgb(ramp(seq(0, 1, length = 1000)), maxColorValue = 255)
  
  N = nrow(conf.mat)
  
  tickdist = ifelse(N > 70, 20, ifelse(N > 30, 10, 5))
  
  if(labels == FALSE){
    low = N - floor(N/tickdist)*tickdist + 1
    x.values = rev(seq(tickdist,N,tickdist))
    y.values = seq(low, N, tickdist)
    lbls = rev(seq(tickdist, N, tickdist))
  }
  else{
    lbls = rownames(conf.mat.ord)
    x.values = 1:N
    y.values = rev(1:N)
  }
  
  lattice::levelplot(t(conf.mat.ord)[,ncol(conf.mat.ord):1], col.regions = colors, 
            xlab = "Loser", ylab = "Winner",
            scales=list(
              x=list(labels=lbls, at = x.values, rot = ifelse(labels == TRUE, 90, 0)),
              y=list(labels=lbls, at = y.values)
            )
  )
}


