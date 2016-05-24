#' Statistical description
#'
#' Provides description of a vector, matrix, data.frame.
#' @param x A data frame, matrix, vector, or formula.
#' @param \dots Additional arguments passed to \code{describe.default}.
#'
#' @examples
#' \dontrun{
#'  describe(turnout)
#'  desc <- describe(turnout)
#'  desc$v1   # print description for just v1
#'  desc[c('v2','v3')]    # print description for two variables.
#'  desc[sort(names(desc))] # print in alphabetic order by column names.
#'
#' # Describing part of a data frame:
#'  with(turnout, describe(v1 ~ v2*v3 + v4) )
#'  with(turnout, describe(~ v2 + v3) )
#'  with(turnout, describe(~ v2 + v3, weights=freqs)) # weighted analysis
#' }
#'
#'
#' @export
`describe` <- function(x, ...) UseMethod("describe")
