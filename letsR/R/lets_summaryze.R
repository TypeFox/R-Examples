#' Summarize variable(s) values in a presence-absence matrix within species' ranges
#' 
#' @author Bruno Vilela & Fabricio Villalobos
#' 
#' @description Based on a Presence-Absence matrix with added variables 
#' (see \code{\link{lets.addvar}}), this function summarizes the values of 
#' such variable(s) per species (across the species' occupied cells. i.e. 
#' within their ranges). 
#' 
#' @param x Presence-absence matrix with variables added.
#' @param pos Column position of the variables of interest.
#' @param xy Logical, if \code{TRUE} the input matrix contains geographic 
#' coordinates in the first two columns. 
#' @param fun Function to be used to summarize the variable per species.
#' 
#' @references Villalobos, F. and Arita, H.T. 2010. The diversity field of 
#' New World leaf-nosed bats (Phyllostomidae). 
#' Global Ecology and Biogeography. 19, 200-211.
#'  
#' @examples \dontrun{
#' data(PAM)
#' data(temp)
#' pamvar <- lets.addvar(PAM, temp)
#' resu <- lets.summarizer(x = pamvar, pos = ncol(pamvar),
#'                         xy = TRUE)
#' }
#' 
#' @seealso \code{\link{lets.addvar}}
#' @seealso \code{\link{lets.addpoly}}
#' @seealso \code{\link{lets.field}}
#' 
#' 
#' @export


lets.summarizer <- function(x, pos, xy = TRUE, fun = mean) {
  
  var <- x[, pos, drop = FALSE]
  sp <- x[, -pos, drop = FALSE]
  
  if (xy) {
    sp <- sp[, -(1:2), drop = FALSE]
  }
  
  Species <- colnames(sp)
  n <- length(Species)
  lpos <- length(pos)
  resum <- matrix(NA, nrow = n, ncol = lpos)
  colnames(resum) <- colnames(var)
  
  for(i in 1:n) {
    vari <- var[(sp[, i] == 1), , drop = FALSE]
    is_all_na <- apply(vari, 2, function(x) {all(is.na(x))})
    if (nrow(vari) == 0 | all(is_all_na)) {
      resum[i, ] <- rep(NA, lpos)
    } else {
      if (any(is_all_na)) {
        resum[i, is_all_na] <- NA
      }
      resum[i, !is_all_na] <- apply(vari[, !is_all_na, drop = FALSE],
                                    2, fun, na.rm = TRUE)
    }
  }
  
  resul <- as.data.frame(cbind(Species, resum))
  
  # Transform in numeric
  for(i in 2:ncol(resul)) {
    resul[, i] <- as.numeric(levels(resul[, i]))[resul[, i]]
  }
  
  return(resul)
}
