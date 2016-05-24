#' Frequency distribution of a variable within a species' range
#' 
#' @author Bruno Vilela
#' 
#' @description Based on a species Presence-Absence matrix including
#'  variables of interest (see \code{\link{lets.addvar}}), the 
#'  function divides a continuous variable into classes and counts 
#'  the frequency of each class within each species' range. 
#' 
#' @param x Presence-absence \code{matrix} with a single variable 
#' added (see \code{\link{lets.addvar}}).
#' @param pos Column number containing the variable of interest.
#' @param groups The number of classes into which the variable will be divided. 
#' Default calculates the number of classes as the 
#' default for a histogram (\code{\link{hist}}). 
#' @param xy Logical, if \code{TRUE} the input matrix contains the geographic 
#' coordinates in the first two columns.
#' 
#' @return A \code{matrix} with species in the rows and the variable's 
#' classes in the columns.
#' 
#' @references Morales-Castilla et al. 2013. Range size patterns of New 
#' World oscine passerines (Aves): insights from differences among migratory
#'  and sedentary clades. Journal of Biogeography, 40, 2261-2273. 
#'
#' @examples \dontrun{
#' data(PAM)
#' data(temp)
#' pamvar <- lets.addvar(PAM, temp)
#' resu <- lets.classvar(x = pamvar, pos = ncol(pamvar), xy = TRUE)
#' }
#'   
#' @export

lets.classvar <- function(x, groups = "default", pos, xy) {
  
  if (xy) {
    sps <- x[, -c(1, 2, pos), drop = FALSE]
  } else {
    sps <- x[, -pos, drop = FALSE]
  }
  
  ric <- ncol(sps)
  
  if (groups == "default") {
    freq <- hist(x[, pos], plot = FALSE)$breaks
    groups <- length(freq) - 1    
  } else {
    freq <- quantile(x[, pos], seq(0, 1, (1 / groups)),
                     na.rm = TRUE)
  }
  
  freqi <- matrix(ncol = groups, nrow = ric)
  
  for(i in 1:ric) {
    freqi[i, ] <- hist(x[(sps[, i] == 1), pos, drop = FALSE],
                       breaks = freq, plot = FALSE)$counts
  }
  
  nomes <- numeric(groups)
  freq <- round(freq, 2)
  
  for(j in 1:groups) {
    nomes[j] <- paste(freq[j], ":", freq[j + 1], sep = "")
  }
  
  rownames(freqi) <- colnames(sps)
  colnames(freqi) <- nomes
  return(freqi)
}
