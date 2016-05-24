#' @title Pop names out
#'
#' @description That is, drop them
#'
#' @export
#'
#' @param .data Input, object of class taxon
#' @param ... Further unnamed args, see examples
#' @return For \code{taxon} inputs gives back a \code{taxon} object. For \code{taxa} inputs
#' gives back a \code{taxa} object. For \code{taxondf} inputs, gives back a \code{taxondf}
#' object.
#' @examples
#' # operating on `taxon` objects
#' out <- make_taxon(genus="Poa", epithet="annua", authority="L.",
#'    family='Poaceae', clazz='Poales', kingdom='Plantae', variety='annua')
#' ## single taxonomic group
#' out %>% pop(family)
#' out %>% pop(genus)
#' out %>% pop(species)
#' ## many taxonomic groups
#' out %>% pop(family, genus, species)
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
#' ## pop out a single taxonomic group
#' df2 %>% pop(order)
#' df2 %>% pop(family)
#' df2 %>% pop(genus)
#'
#' ## pop out many taxonomic groups
#' df2 %>% pop(order, family)
#' df2 %>% pop(order, genus)
#'
#' # From taxa object
#' df2 %>% scatter %>% pop(family)
#' df2 %>% scatter %>% pop(family, species)
#' df2 %>% scatter %>% pop(family, species, genus)
pop <- function(.data, ...) {
  UseMethod("pop")
}

#' @export
pop.taxon <- function(.data, ...){
  tmp <- .data$grouping
  name <- vars(...)
  taxon(binomial = .data$binomial,
        grouping = do.call("grouping", tmp[!names(tmp) %in% name]))
}

#' @export
pop.taxa <- function(.data, ...){
  taxa(lapply(.data, pop, ...))
}

#' @export
pop.taxondf <- function(.data, ...){
  var <- vars(...)
  check_vars(var, names(.data))
  select_(.data, .dots = paste("-", var, sep = ""))
}
