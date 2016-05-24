#' Sort taxon or taxondf objects by one or more choices
#'
#' @import lazyeval
#' @export
#' @param x Input, object of class taxondf
#' @param ... Comma separated list of unquoted expressions. You can treat variable names
#' like they are positions. Use positive values to select variables; use negative values
#' to drop variables.
#' @param .dots Use sort_() to do standard evaluation
#' @examples
#' # operating on taxonomic data.frames
#' df <- data.frame(class=c('Magnoliopsida','Magnoliopsida','Magnoliopsida',
#'                          'Magnoliopsida','Magnoliopsida','Magnoliopsida'),
#'          order=c('Asterales','Asterales','Fagales','Poales','Poales','Poales'),
#'          family=c('Asteraceae','Asteraceae','Fagaceae','Poaceae','Poaceae','Poaceae'),
#'          genus=c('Helianthus','Helianthus','Quercus','Poa','Festuca','Holodiscus'),
#'          stringsAsFactors = FALSE)
#' (df2 <- taxon_df(df))
#'
#' ## sort the taxonomic data.frame
#' df2 %>% sort(order)
#' df2 %>% sort(family)
#' df2 %>% sort(genus)
#' df2 %>% sort(genus, order)
#'
#' ## reverse order
#' df2 %>% sort(desc(genus))
sort <- function(x, ...) {
  sort_(x, .dots = lazyeval::lazy_dots(...))
}

#' @export
#' @rdname sort
sort_ <- function(x, ..., .dots) {
  UseMethod("sort_")
}

#' @export
#' @rdname sort
sort_.taxondf <- function(x, ..., .dots) {
  dots <- lazyeval::all_dots(.dots, ..., all_named = TRUE)
  structure(as.data.frame(arrange_impl(x, dots)), class = c('taxondf', 'data.frame'))
}

arrange_impl <- function(data, dots) {
  .Call("dplyr_arrange_impl", PACKAGE = "dplyr", data, dots)
}

#' @export
#' @rdname sort
desc <- function(x) {
  -xtfrm(x)
}
