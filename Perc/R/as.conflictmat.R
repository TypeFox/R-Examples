# transform an edgelist into a matrix
# 
# @param edgelist a 2-column (or 3-column for weighted edgelist) dataframe/matrix of edges. The winner is in the 1st column by default. For weighted edgelist, the third column should be the weight. 
# @param weighted If the edgelist is a 3-column weighted edgelist, use \code{weighted = TRUE}. 
# @param swap.order If the winner is in the 2nd column, specify as \code{TRUE}.
# @return a named matrix with \code{[i,j]}th entry equal to the number of times \code{i} wins over \code{j}. 
# It is the matrix representation of the edgelist.
# 
# @seealso \code{\link{conductance}}

edgelisttomatrix <- function(edgelist, weighted = FALSE, swap.order = FALSE) {
  
  if (ncol(edgelist) > 3) {
    stop("edgelist should be of 2 column, or 3-column for weighted edgelist")
  }
  
  if (any(lapply(edgelist[, 1:2], class) == 'factor')){
    edgelist[, 1:2] <- apply(edgelist[, 1:2], 2, as.character)
    warning("Factor vector(s) in your data into were automatically converted into character vector(s).")
  }
  
  if (swap.order == TRUE){
    edgelist[, 1:2] <- edgelist[, 2:1]
  }
  
  if (any(edgelist[,1] == edgelist[,2])) {
    rowIndex <- which(edgelist[,1] == edgelist[,2])
    edgelist <- edgelist[-rowIndex, ]
    warning(
      paste("Row number", 
            paste(rowIndex, collapse = ","), 
            "in your raw data is removed because the initiator and recipient are the same. Self-look is not allowed for this method.")
      )
  }
  
  subjects = unique(sort(as.matrix(edgelist[,1:2]))) # work better for IDs of character
  # subjects = sort(unique(c(edgelist[,1], edgelist[,2])))
  N = length(subjects)
  if (N > 10000){
    stop("No more than 10000 unique subjects.")
  }
  
  mat = matrix(0, N, N)
  
  
  if (weighted == TRUE){
    
    if (ncol(edgelist) != 3){
      stop("Input a matrix or dataframe with three columns, with the third column being Frequency of the interaction")
    }
    
    if (anyDuplicated(edgelist[,1:2]) != 0) {
      warning(
        "dyads in the weighted edgelist are not unique; the sum of frequencies is taken for duplicated rows."
        )
      edgelist <- sumDuplicate(edgelist)
    }
    
    
    # transform the weighted edgelist into a matirx
    
    for(i in 1:nrow(edgelist)){
      subject1 = which(subjects == edgelist[i,1])
      subject2 = which(subjects == edgelist[i,2])
      mat[subject1, subject2] = edgelist[i, 3]
    }
    
  } else {
    
    if (ncol(edgelist) != 2){
      stop("edgelist should be a dataframe or matrix of two columns. If it is a weighted edgelist, it should be a matrix or dataframe of 3 columns and use the argument 'weighted = TRUE'")
    }
    
    for(i in 1:nrow(edgelist)){
      subject1 = which(subjects == edgelist[i,1])
      subject2 = which(subjects == edgelist[i,2])
      mat[subject1, subject2] = mat[subject1, subject2] + 1
    }
  }
  
  rownames(mat) = subjects
  colnames(mat) = subjects
  
  return(mat)
}


#' convert to a matrix of \code{conf.mat} class
#' 
#' \code{as.conflictmat} convert an edgelist or a win-loss raw matrix to a matrix of \code{conf.mat} class
#' @param Data either a dataframe or a matrix, representing raw win-loss interactions using either an edgelist or a matrix. 
#' By default, winners are represented by IDs in the 1st column for an edgelist, and by row IDs for a matrix. 
#' Frequency of interactions for each dyad can be represented either by multiple occurrences of the dyad for a 2-column edgelist, or
#' by a third column specifying the frequency of the interaction for a 3-column edgelist.
#' @param swap.order If the winner is placed in the 2nd column for an edgelist or as the column name for a matrix, specify as \code{TRUE}. By default, winners are placed in the first column of an edgelist or in the row names of a matrix.
#' @param weighted If the edgelist is a 3-column edgelist in which weight was specified by frequency, use \code{weighted = TRUE}. 
#' @return a named matrix with the \code{[i,j]}th entry equal to the number of times \code{i} wins over \code{j}.
#' @details \code{conf.mat} is short for "Conflict Matrix". \code{conf.mat} is 
#' a class of R objects. It is required to use \code{as.conflictmat} to convert your
#' raw edgelist or raw win-loss matrix into a matrix of \code{conf.mat} object before
#' using other functions to find (in)direct pathways and computing dominance probabilities.
#' 
#' Note, when using a 3-column edgelist (e.g. a weighted edgelist) to represent raw win-loss interactions, each dyad must be unique. If more than one rows are found with the same initiator and recipient,
#' sum of the frequencies will be taken to represent the freqency of interactions between this unique dyad. A warning message will prompt your attention to the accuracy of your raw data when duplicate dyads were found in a three-column edgelist.
#' 
#' 
#' @seealso \code{\link{findIDpaths}}, \code{\link{countPaths}}, \code{\link{transitivity}}, \code{\link{conductance}}
#' 
#' @examples
#' confmatrix <- as.conflictmat(sampleEdgelist, swap.order = FALSE)
#' confmatrix2 <- as.conflictmat(sampleRawMatrix, swap.order = FALSE)
#' confmatrix3 <- as.conflictmat(sampleWeightedEdgelist, weighted = TRUE, swap.order = FALSE)
#' @export

as.conflictmat = function(Data, weighted = FALSE, swap.order = FALSE){
  if (ncol(Data) > 3 & ncol(Data) != nrow(Data)) {
    stop("check your raw data: A edgelist should be of either 2 or 3 columns. If it is a win-loss matrix, the column number should be equal to row number.")
  }
  
  if (ncol(Data) == nrow(Data)){
    # if values on diagonal are not all zeros, return warnings.
    if (any(diag(as.matrix(Data)) != 0)){
      index <- which(diag(as.matrix(Data)) != 0)
      mat <- as.matrix(Data)
      diag(mat)[index] <- 0
      warning(paste("Non-zero values in your raw win-loss matrix at Row", 
                    paste(index, collapse = ","), 
                    "and column", 
                    paste(index, collapse = ","), 
                    " are converted to 0; Non-zero values are not allowed on the diagonal in your raw win-loss matrix."))
    }
    if (swap.order == TRUE) {
      mat <- t(as.matrix(Data))
    } else{
      mat <- as.matrix(Data)
      
    }
    # update 2016.1.18: sorted matrix by colnames
    sorted_subjects = unique(sort(colnames(mat)))
    mat_sorted = mat[sorted_subjects, sorted_subjects]
  } else {
    mat_sorted <- edgelisttomatrix(Data, weighted, swap.order)
  }
  class(mat_sorted) = c("conf.mat", "matrix")
  return(mat_sorted) 
}


#### internal functions


sumDuplicate <- function(weightedEdgelist) {
  uniqueEdgelist <- unique(weightedEdgelist[,1:2])
  for (i in 1:nrow(uniqueEdgelist)){
    uniqueEdgelist[i,3] <- 
      sum(
        weightedEdgelist[
          match.2coldf(weightedEdgelist[,1:2],  uniqueEdgelist[i,]),
          3])
  }
  names(uniqueEdgelist) <- names(weightedEdgelist)
  return(uniqueEdgelist)
}


match.2coldf <- function(dataframe, values) {
  # dataframe should be of two columns
  # values should be a vector of length 2, or a row of dataframe of two columns
  rowIndex <- intersect(which(dataframe[,1] == values[[1]]),
                        which(dataframe[,2] == values[[2]]))
  return(rowIndex)
}