#'A cleaning function for special characters
#'
#' This function is a helper for cleaning xtable outputs in preperation for PDF production
#'
#' @param str The text to be sanitised
#' 
#' 
#' @keywords knitr pdflatex generate PDF Rnw
#' @family helper  
#' 
#' @examples
#' sanitise('[&%#<>\\')
#' 
#' @export



sanitise <- function(str) {
    result <- str
    result <- gsub("\\\\", "SANITIZE.BACKSLASH", result)
    result <- gsub("$", "\\$", result, fixed = TRUE)
    result <- gsub(">", "$>$", result, fixed = TRUE)
    result <- gsub("<", "$<$", result, fixed = TRUE)
    result <- gsub("|", "$|$", result, fixed = TRUE)
    result <- gsub("{", "\\{", result, fixed = TRUE)
    result <- gsub("}", "\\}", result, fixed = TRUE)
    result <- gsub("[", "{}[", result, fixed = TRUE)  # Line added, everything else is taken from xtable function definition
    result <- gsub("%", "\\%", result, fixed = TRUE)
    result <- gsub("&", "\\&", result, fixed = TRUE)
    result <- gsub("_", "\\_", result, fixed = TRUE)
    result <- gsub("#", "\\#", result, fixed = TRUE)
    result <- gsub("\u00a3", "\\pounds ", result, fixed = TRUE)
    result <- gsub("^", "\\verb|^|", result, fixed = TRUE)
    result <- gsub("~", "\\~{}", result, fixed = TRUE)
    result <- gsub("SANITIZE.BACKSLASH", "$\\backslash$", result, fixed = TRUE)
    return(result)
} 
