#' Parse taxon or taxondf objects by a range of names
#'
#' @export
#' @param .data Input, object of class taxon
#' @param ... Pass in two unquoted taxonomic rank names, and only two. May make this
#' more flexible in the future.
#' @return A single or list of \code{taxon} class objects
#' @examples
#' # operating on `taxon` objects
#' out <- make_taxon(genus="Poa", epithet="annua", authority="L.",
#'    family='Poaceae', clazz='Poales', kingdom='Plantae', variety='annua')
#' out %>% span(kingdom, genus)
#'
#' # operating on taxonomic data.frames
#' df <- data.frame(class=c('Magnoliopsida','Magnoliopsida','Magnoliopsida',
#'                          'Magnoliopsida','Magnoliopsida','Magnoliopsida'),
#'          order=c('Asterales','Asterales','Fagales','Poales','Poales','Poales'),
#'          family=c('Asteraceae','Asteraceae','Fagaceae','Poaceae','Poaceae','Poaceae'),
#'          genus=c('Helianthus','Helianthus','Quercus','Poa','Festuca','Holodiscus'),
#'          stringsAsFactors = FALSE)
#' (df2 <- taxon_df(df))
#'
#' ## filter to get a range of classes
#' df2 %>% span(order, genus)
#' df2 %>% span(family, genus)
#'
#' ## from taxa object
#' df2 %>% scatter %>% span(family, species)
span <- function(.data, ...) {
  UseMethod("span")
}

#' @export
span.taxon <- function(.data, ...) {
  var <- vars(...)
  taxonparse(.data, var)
}

#' @export
span.taxa <- function(.data, ...) {
  lapply(.data, span, ...)
}

#' @export
span.taxondf <- function(.data, ...) {
  var <- vars(...)
  if(length(var) > 2) stop("Pass in only two rank names", call. = FALSE)
  check_vars(var, names(.data))
  matches <- sapply(var, grep, x=names(.data))
  .data[fill_nums(matches)]
}

taxonparse <- function(w, vars){
  tmp <- w$grouping
  if(length(vars) != 2) stop("Pass in only two rank names", call. = FALSE)
  check_vars(vars, names(tmp))
  matches <- sapply(vars, grep, x=names(tmp))
  w$grouping <- do.call("grouping", tmp[fill_nums(matches)])
  return(w)
}
