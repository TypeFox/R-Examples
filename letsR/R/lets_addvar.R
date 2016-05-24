#' Add variables (in raster format) to a PresenceAbscence object
#' 
#' @author Bruno Vilela
#' 
#' @description Add variables (in raster format), usually environmental, to a PresenceAbsence object. Variables are included as additional columns containing the aggregate/summarize value of the variable(s) in each cell of the presence-absence matrix.
#'
#' @param x A \code{\link{PresenceAbsence}} object. 
#' @param y Variables to be added in \code{Raster} or \code{RasterStack} format.
#' @param onlyvar If \code{TRUE} only the \code{matrix} object will be returned.
#' @param fun Function used to aggregate the variables(s) values over each cell.
#' Note that this will only work for variables with a resolution value smaller
#' (i.e. higher resolution) than the PAM.
#' 
#' @note The \code{PresenceAbsence} and the \code{Raster} 
#' variable must be in the same projection.
#' 
#' @return The result is a presence-absence matrix of species with 
#' the variables added as columns at the right-end of the matrix 
#' (but see the 'onlyvar' argument).
#'  
#' @seealso \code{\link{lets.presab.birds}}
#' @seealso \code{\link{lets.presab}}
#' @seealso \code{\link{lets.addpoly}}
#' 
#' @examples \dontrun{
#' data(temp)  # Global mean temperature
#' data(PAM)  # Phyllomedusa presence-absence matrix
#' # Mean temperature
#' PAM_temp_mean <- lets.addvar(PAM, temp)
#' # Standard deviation of temperature
#' PAM_temp_sd <- lets.addvar(PAM, temp, fun = sd, onlyvar = TRUE)
#' # Mean and SD in the PAM
#' PAM_temp_mean_sd <- cbind(PAM_temp_mean, PAM_temp_sd)
#' }
#' 
#' @export

lets.addvar <- function(x, y, onlyvar = FALSE, fun = mean) {
  
  # Error control for different projections
  if (projection(x[[2]]) != projection(y)) {
    stop("The PAM object and the variable (raster) 
         should be in the same projection")
  }
  
  # Crop variable by the PAM
  var_c <- crop(y, x[[2]])
  
  # Check the resolution, who has the bigger
  res1 <- res(var_c)[1]
  res2 <- res(x[[2]])[1]
  if (res2 > res1) {
    var_a <- aggregate(var_c, fact = res2 / res1,
                       na.rm = TRUE, fun)
  }
  if (res2 < res1) {
    var_a <- disaggregate(var_c, fact = res1 / res2,
                          na.rm = TRUE)
  }
  if (res2 == res1) {
    var_a <- var_c
  }
  
  # Adjust rasters
  var_r <- resample(var_a, x[[2]])
  
  # Extract values
  var_e <- extract(var_r, x[[1]][, 1:2]) 
  var_e <- as.matrix(var_e)
  
  # Set the names
  colnames(var_e) <- paste(names(y), as.character(substitute(fun)),
                           sep = "_")
  
  # Return result
  resultado <- cbind(x[[1]], var_e)
  if (onlyvar) {
    return(var_e) 
  } else {
    return(resultado)
  }
}
