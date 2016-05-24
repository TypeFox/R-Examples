#' Subset a PresenceAbsence object based on species names
#' 
#' @author Bruno Vilela
#' 
#' @description Subset a PresenceAbsence object based on species character vector 
#' provided by the user.
#' 
#' @param x A \code{\link{PresenceAbsence}} object.
#' @param names Character vector with species names to subset the \code{PresenceAbsence} object.
#' @param remove.cells Logical, if \code{TRUE} the final matrix will not contain cells in the 
#' grid with a value of zero (i.e. sites with no species present).
#' 
#' 
#' @return The result is an object of class PresenceAbsence subseted.
#' 
#' 
#' @seealso \code{\link{plot.PresenceAbsence}}
#' @seealso \code{\link{lets.presab.birds}} 
#' 
#' @examples \dontrun{
#' data(PAM)
#' # PAM before subset
#' plot(PAM, xlab = "Longitude", ylab = "Latitude",
#'      main = "Phyllomedusa species richness")
#' 
#' # Subset PAM to the first 20 species
#' PAMsub <- lets.subsetPAM(PAM, PAM[[3]][1:20])
#' plot(PAMsub, xlab = "Longitude", ylab = "Latitude",
#'      main = "Phyllomedusa species richness")
#' }
#' 
#' @export

lets.subsetPAM <- function(x, names, remove.cells = TRUE) {
  
  if (class(x) != "PresenceAbsence") {
    stop("x argument must be a PresenceAbsence object")
  }
  
  if (class(names) != "character") {
    stop("names argument must be a character object")
  }
  
  if (class(remove.cells) != "logical") {
    stop("remove.cells argument must be a TRUE or FALSE")
  }
  
  
  # Get species position name
  pos <- colnames(x[[1]]) %in% names  
  
  errorcont <- names %in% colnames(x[[1]]) 
  
  if (!any(errorcont)) {
    stop("None of the names provided match with PAM species")
  }
  
  if (any(!errorcont)) {
    warning(paste("One or more names",
                  "provided, do not",
                  "match any of the",
                  "PAM species"))
  }
  
  x[[1]] <- x[[1]][, c(1:2, which(pos)), drop = FALSE]
  
  if (remove.cells) {
    x[[1]] <- .removeCells(x[[1]])
  }
  
  rich <- rowSums(x[[1]][, -(1:2), drop = FALSE])
  x[[2]] <- rasterize(x[[1]][, c(1:2), drop = FALSE], x[[2]], rich)
  x[[3]] <- x[[3]][pos]
  
  return(x)
}