
#' Convert municipality names into standard versions harmonized within the package
#'
#' @param municipality.names municipality names to convert
#' @return standardized municipality names
#' @references
#' See citation("sorvi") 
#' @author Leo Lahti \email{louhos@@googlegroups.com}
#' @examples \dontrun{tmp <- convert_municipality_names("Koski.Tl")}
#' @export
#' @keywords internal
convert_municipality_names <- function (municipality.names) {

  f <- system.file("extdata/municipality_synonymes.csv", package = "sorvi")
  syn <- read.csv(f, sep = "\t")		 
  
  harmonized.names <- c()
  for (nam in municipality.names) {
    if (nam %in% syn$synonyme) {
      harmonized.names[[nam]] <- syn[which(syn$synonyme == nam), "name"]
    } else {
      harmonized.names[[nam]] <- nam      
    }
  }
  
  harmonized.names

}

