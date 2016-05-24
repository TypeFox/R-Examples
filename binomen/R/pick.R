#' Pick names
#'
#' @export
#' @param .data Input, object of class taxon
#' @param ... Further unnamed args, see examples
#' @return For \code{taxon} inputs, gives back a \code{taxonref} object. For \code{taxondf}
#' inputs, gives back \code{taxondf}.
#' @examples
#' # operating on `taxon` objects
#' out <- make_taxon(genus="Poa", epithet="annua", authority="L.",
#'    family='Poaceae', clazz='Poales', kingdom='Plantae', variety='annua')
#' out %>% pick(family)
#' out %>% pick(genus)
#' out %>% pick(species, genus)
#' out %>% pick(species) %>% name()
#' out %>% pick(species) %>% uri()
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
#' ## select single or many taxonomic classes
#' df2 %>% pick(order)
#' df2 %>% pick(family, genus)
#'
#' # From taxa object
#' df2 %>% scatter %>% pick(family)
#' df2 %>% scatter %>% pick(family, species)
#' df2 %>% scatter %>% pick(family, species, genus)
pick <- function(.data, ...) {
  UseMethod("pick")
}

#' @export
pick.taxon <- function(.data, ...){
  tmp <- .data$grouping
  name <- vars(...)
  taxon(binomial = .data$binomial,
        grouping = do.call("grouping", tmp[names(tmp) %in% name]))
}

#' @export
pick.taxa <- function(.data, ...){
  taxa(lapply(.data, pick, ...))
}

#' @export
pick.taxondf <- function(.data, ...){
  var <- vars(...)
  check_vars(var, names(.data))
  select_(.data, .dots = var)
}

# helpers ---------------------------
fill_nums <- function(x) seq(from=min(x), to=max(x), by=1)

vars <- function(...) as.character(dots(...))

dots <- function(...){
  eval(substitute(alist(...)))
}

check_vars <- function(x, y){
  if( !all(x %in% y) ) stop(sprintf("%s not a valid taxonomic rank in your dataset", paste0(x[!x %in% y], collapse = ", ")), call. = FALSE)
}
