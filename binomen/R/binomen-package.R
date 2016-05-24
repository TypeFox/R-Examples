#' Taxonomic class specification and parsing methods
#'
#' @importFrom methods is
#' @importFrom stats setNames
#' @import lazyeval
#' @importFrom dplyr select_ rbind_all
#' @importFrom jsonlite toJSON
#' @name binomen-package
#' @aliases binomen
#' @docType package
#' @keywords package
#' @author Scott Chamberlain \email{myrmecocystus@@gmail.com}
#'
#' @examples
#' library("binomen")
#'
#' # operating on `taxon` objects
#' out <- make_taxon(genus="Poa", epithet="annua", authority="L.",
#'    family='Poaceae', clazz='Poales', kingdom='Plantae', variety='annua')
#' # get single name
#' out %>% pick(family)
#' out %>% pick(genus)
#' out %>% pick(species)
#' out %>% pick(species) %>% name()
#' out %>% pick(species) %>% uri()
#' # get range of names
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
#' ## select single taxonomic class
#' df2 %>% pick(order)
#' df2 %>% pick(family, genus)
#'
#' ## filter to get a range of classes
#' df2 %>% span(order, genus)
#' df2 %>% span(family, genus)
NULL

#' Lookup-table for IDs of taxonomic ranks
#'
#' @name rank_table
#' @docType data
#' @keywords data
#' @return data.frame with two columns and 34 rows. The two columns:
#' \itemize{
#'  \item rankid integer; smaller numbers are higher ranks.
#'  \item ranks character; rank name. Some rows have more than one name,
#'  in which the names are considered of equal height.
#' }
NULL
