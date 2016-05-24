
#' Conversions between municipality codes and names
#'
#' @param ids NULL 
#' @param municipalities NULL 
#'
#' @return Depending on the input. Converted id or name vector, or full conversion table.
#' @export 
#' @references
#' See citation("sorvi") 
#' @author Leo Lahti \email{louhos@@googlegroups.com}
#' @examples  \dontrun{conversion.table <- convert_municipality_codes()}
#' @keywords utilities

convert_municipality_codes <- function (ids = NULL, municipalities = NULL) {
 
  # Reading municipality information from the web
  df <- get_municipality_info_mml()	

  conversion.table <- df[, c("Kunta", "Kunta.FI")]
  names(conversion.table) <- c("id", "name")

  #write.csv(conversion.table, file = "../inst/extdata/conversiontable.tab", quote = FALSE, row.names =FALSE)
  #conversion.table <- read.csv(paste(system.file("extdata", package = "sorvi"), 
  # 		     	"/conversiontable.tab", sep = ""))

  conversion.table$id <- as.character(conversion.table$id)
  conversion.table$name <- as.character(conversion.table$name)

  res <- conversion.table

  if (!is.null(ids)) {
    res <- conversion.table$name[match(as.character(ids), conversion.table$id)]
    names(res) <- ids
  } else if (!is.null(municipalities)) {
    res <- conversion.table$id[match(as.character(municipalities), conversion.table$name)]
    names(res) <- municipalities
  } 

  res

}


