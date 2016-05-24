#' Scatter each taxon in a taxondf to a taxon object
#'
#' @export
#'
#' @param x A taxonomic data.frame
#' @param ... Further args, ignored for now
#' @return Gives a \code{taxa} object, with each individual component a row from your
#' data.frame, and of class \code{taxon}
#' @details Right now, \code{assemble} may not give back the identical data.frame that one
#' would pass to \code{scatter}.
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
#' ## scatter each taxon into a taxon class
#' df2 %>% scatter()
#'
#' ## re-assemble
#' df2
#' df2 %>% scatter()
#' df2 %>% scatter() %>% assemble
scatter <- function(x, ...) {
  UseMethod("scatter")
}

#' @export
scatter.taxondf <- function(x, ...) {
  x <- class2clazz(x)
  taxa(unname(apply(x, 1, function(y) {
    do.call("make_taxon", as.list(y))
  })))
}

class2clazz <- function(x){
  if ("class" %in% names(x)) {
    names(x)[which(names(x) == "class")] <- "clazz"
    x
  } else {
    x
  }
}

#' @export
#' @rdname scatter
assemble <- function(x, ...) {
  UseMethod("assemble")
}

#' @export
#' @rdname scatter
assemble.taxa <- function(x, ...) {
  tmp <- lapply(x, "[[", "grouping")
  x <- as.data.frame(rbind_all(lapply(tmp, function(b) {
    data.frame(lapply(b, function(n) setNames(n[['name']], n[['rank']])),
               stringsAsFactors = FALSE, row.names = NULL)
  })))
  taxon_df(x)
}
