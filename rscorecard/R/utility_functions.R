#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @importFrom jsonlite fromJSON
#' @importFrom dplyr filter_
#' @export
dplyr::filter_
#' @importFrom lazyeval interp
#' @export
lazyeval::interp

## paste pipes
`%+%`  <- function(a,b) paste(a, b, sep = '')
`%+|%` <- function(a,b) paste(a, b, sep = '|')
`%+&%` <- function(a,b) paste(a, b, sep = '&')

