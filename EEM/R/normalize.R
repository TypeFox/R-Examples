#' Normalize data
#' 
#' Normalize data (area under the curve = 1)
#' 
#' @param EEM_uf Unfolded EEM matrix where columns are wavelength condition and rows are samples
#' 
#' @return A matrix of normalized data
#' 
#' @details The unfolded EEM data can be normalized by dividing each variable by the 
#' sum of the absolute value of all variables in a sample, such that the summation 
#' of absolute values of all variables in each sample was equal to 1. This is 
#' can be used to reduce the scaling difference,which is common in spectroscopic 
#' applications. This difference is usually caused by the scattering effect, 
#' source/detector variation and instrumental sensitivity.
#' 
#' @examples
#' data(applejuice)
#' applejuice_uf <- unfold(applejuice) # unfold list into matrix
#' applejuice_uf_norm <- normalize(applejuice_uf) # normalize data
#' 
#' rowSums(abs(applejuice_uf_norm), na.rm = TRUE) # the absolute sum of each row equal to 1
#' 
#'   
#' @export
#' 
normalize <-
function(EEM_uf){
  numVar <- dim(EEM_uf)[2] # number of variables
  
  # calculate norm where norm = sum(abs(data))
  normTmp <- as.matrix(rowSums(abs(EEM_uf), na.rm = TRUE))
  normTmpMatrix <- do.call("cbind", rep(list(normTmp), numVar))
  
  # find norm =0 if there is any norm = 0 then stop function
  stopifnot(sum(normTmp == 0) == 0)
  
  # normalize data
  EEM_uf_norm <- EEM_uf / normTmpMatrix
  return(EEM_uf_norm)
}
