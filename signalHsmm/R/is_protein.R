#' Protein test
#'
#' Checks if an object is a protein (contains letters from one-letter amino acid code).
#'
#' @param object \code{character} vector where each elemenents represent one amino acid.
#' @return \code{TRUE} or \code{FALSE}.
#' @export

is_protein <- function(object) {
  #only amino acids
  only_aa <- all(toupper(object) %in% c(a()[-1], "X", "J", "Z", "B", "U"))
  #only nucleotides
  only_nuc <- all(toupper(object) %in% c("A", "C", "G", "T"))
 
  only_aa & !only_nuc
}
