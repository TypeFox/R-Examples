#' @title Common regular expression patterns used by lettercase 
#' 
#' @description These are patterns that are used in the lettercase package.
#' 
#' @details 
#'   \code{lettercase} provides a number of regular expression patterns that are
#'   used by the case conversion functions.
#'   
#'   These are not exported; to use them prefix them with \code{lettercase:::*}
#'   
#' @rdname patterns
  pattern_whitespace <- '\\s'
 
#' @rdname patterns 
  pattern_whitespace_like <- "[-_]"

#' @rdname patterns
  pattern_separators <- '[\\s_\\.-]'

#' @rdname patterns 
  pattern_ucfirst <- '\\b([a-z])'

#' @rdname patterns
  pattern_word       <- '\\w'

#' @rdname patterns
  pattern_nonword    <- '\\W'

#' @rdname patterns
  pattern_uppercase <- '[A-Z]'

#' @rdname patterns
  pattern_lowercase <- '[a-z]'
  