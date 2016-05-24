#' @title Pick out parts by name
#'
#' @description This suite of functions act on taxon or taxonref objects,
#' and pick out object elements by the name of the function.
#'
#' @name parts
#'
#' @param x Input, object of class taxon or taxonref
#' @param unname (logical) Unname output elements? Ignored when input is of class
#' \code{taxonref}. Default: \code{TRUE}
#' @return For \code{taxon} inputs, gives back a \code{taxonref} object. For \code{taxondf}
#' inputs, gives back \code{taxondf}.
#' @examples
#' # operating on `taxon` objects
#' out <- make_taxon(genus="Poa", epithet="annua", authority="L.",
#'    family='Poaceae', clazz='Poales', kingdom='Plantae', variety='annua')
#'
#' out %>% name()
#' out %>% uri()
#' out %>% rank()
#' out %>% taxonid()
#'
#' ## or don't unname the output
#' out %>% name(unname = FALSE)
#'
#' # operating on `taxonref` objects
#' res <- taxonref("genus", "Poa", 56, "http://scottchamberlain.info/")
#' res %>% name()
#' res %>% uri()
#' res %>% rank()
#' res %>% taxonid()

#' @export
#' @rdname parts
name <- function(x, unname = TRUE) {
  UseMethod("name")
}

#' @export
#' @rdname parts
name.taxon <- function(x, unname = TRUE) {
  na_me(pluck(x$grouping, "name", ""), unname)
}

#' @export
#' @rdname parts
name.taxonref <- function(x, unname = TRUE) {
  x$name
}


#' @export
#' @rdname parts
uri <- function(x, unname = TRUE) {
  UseMethod("uri")
}

#' @export
#' @rdname parts
uri.taxon <- function(x, unname = TRUE) {
  na_me(pluck(x$grouping, "uri", ""), unname)
}

#' @export
#' @rdname parts
uri.taxonref <- function(x, unname = TRUE) {
  x$uri
}


#' @export
#' @rdname parts
rank <- function(x, unname = TRUE) {
  UseMethod("rank")
}

#' @export
#' @rdname parts
rank.taxon <- function(x, unname = TRUE) {
  na_me(pluck(x$grouping, "rank", ""), unname)
}

#' @export
#' @rdname parts
rank.taxonref <- function(x, unname = TRUE) {
  x$uri
}


#' @export
#' @rdname parts
taxonid <- function(x, unname = TRUE) {
  UseMethod("taxonid")
}

#' @export
#' @rdname parts
taxonid.taxon <- function(x, unname = TRUE) {
  na_me(pluck(x$grouping, "id", ""), unname)
}

#' @export
#' @rdname parts
taxonid.taxonref <- function(x, unname = TRUE) {
  x$id
}


# helper fxns --------------
na_me <- function(x, z) {
  if (z) {
    unname(x)
  } else {
    x
  }
}
