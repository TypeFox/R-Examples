#' Binarize digital PCR data
#' 
#' Transforms multinomial (number of molecules per partition) or continuous (fluorescence)
#' digital PCR data to binary (positive/negative partition) format.
#' 
#' @aliases binarize
#' @param input object of the class \code{\linkS4class{adpcr}} or
#' \code{\linkS4class{ddpcr}} with one of following types:\code{"ct"}, \code{"fluo"} or
#' \code{"nm"}.
#' @return object of the class \code{\linkS4class{adpcr}} or
#' \code{\linkS4class{ddpcr}} (depending on \code{input}) with type \code{"np"}.
#' @author Michal Burdukiewicz.
#' @keywords manip
#' @export
#' @examples
#' 
#' #adpcr object
#' rand_array <- sim_adpcr(200, 300, 100, pos_sums = FALSE, n_panels = 1)
#' binarize(rand_array)
#' 
#' #ddpcr object
#' rand_droplets <- sim_ddpcr(200, 300, 100, pos_sums = FALSE, n_exp = 1)
#' binarize(rand_droplets)
binarize <- function(input) {
  if (class(input) %in% c("adpcr", "ddpcr")) {
    if(slot(input, "type") %in% c("tp", "tnp"))
      stop("Cannot binarize already binary data.")
    positive_threshold <- if (class(input) == c("adpcr")) {
      slot(input, "breaks")[2]
    } else {
      slot(input, "threshold")
    }
  } else {
    stop("Input must have 'adpcr' or 'ddpcr' class.")
  }
  
  bin_data <- slot(input, ".Data") >= positive_threshold
  storage.mode(bin_data) <- "integer"
  slot(input, ".Data") <-  bin_data
  slot(input, "type") <- "np"
  if (class(input) == c("adpcr")) {
    slot(input, "breaks") <- c(0, 1)
  } else {
    slot(input, "threshold") <- 1
  }
  input
}