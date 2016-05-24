#' @title harmonize_names
#' @description Harmonize names
#'
#' @param x A character vector 
#' @param synonymes synonyme table with the fields 'synonyme' and 'name'
#' @param remove.unknown Logical. Remove terms that do not have synonymes.
#' @param check.synonymes Check the synonyme table
#'
#' @return Harmonized vector where synonymes are renamed by the selected names
#'
#' @export
#'
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("sorvi")
#' 
#' @examples \dontrun{x2 <- harmonize_names(x, synonymes)}
#' @keywords utilities
harmonize_names <- function (x, synonymes, remove.unknown = FALSE, check.synonymes = TRUE) {

  x <- as.character(x)

  # Check which terms are not on the synonyme list and add them there		
  if (!remove.unknown) {
    r <- setdiff(x, synonymes$synonyme)
    synonymes <- rbind(synonymes,
    	      as.data.frame(list(name = r, synonyme = r)))    
  }

  # Polish the synonyme table
  if (check.synonymes) {
    synonymes <- check_synonymes(synonymes)
  }
  
  # Map synonymes to selected names: NA if mapping not available
  # FIXME: speed up by considering unique x only, then remapping
  xx <- c()
  for (i in 1:length(x)) {
    xh <- unique(as.character(synonymes$name[synonymes$synonyme == x[[i]]]))
    if (length(xh) == 1) {
      xx[[i]] <- xh
    } else {
      warning(paste("No unique mapping available for", x[[i]]))
      xx[[i]] <- NA
    }
  }
  
  # message("Return data frame")
  data.frame(list(name = xx, original = x))

}


