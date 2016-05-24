#' @include dataset.R
#' @include dataset-characteristics.R
{}



#' Dataset characterization framework
#'
#' Implementation of a map/reduce approach to characterize a dataset
#' with given dataset characteristics.
#'
#' @param x A \code{\link[=as.dataset]{dataset}} object
#' @param y A \code{DatasetCharacteristics} object; e.g.,
#'   \code{\link{StatlogCharacteristics}}
#' @param verbose Show information during execution
#' @param index Characterize only a subset
#' @param ... Ignored
#'
#' @return The characterization matrix (1 row and as many columns as
#'   characteristics
#'
#' @examples
#'   data("iris")
#'   ds <- as.dataset(Species ~ ., iris)
#'   characterize(ds, StatlogCharacteristics)
#'
#' @family dataset-characterization
#'
#' @references
#'   See \emph{Eugster et al. (2010)} in \code{citation("benchmark")}.
#'
#' @export
characterize <- function(x, y, verbose = FALSE, index = NULL, ...) {
  stopifnot(is(x, 'dataset'))
  stopifnot(is(y, 'characteristics'))

  stopifnot(y$requirements())

  d <- map(x, y, verbose = verbose, index = index)
  d <- reduce(d, y, verbose = verbose)

  d <- as.matrix(as.data.frame(d))
  #rownames(d) <- x$dataname()

  d
}



### Characterization map/reduce framework: ###########################

map <- function(x, y, ...) {
  UseMethod('map')
}


map.dataset <- function(x, y, verbose = TRUE, index = NULL, ...) {
  stopifnot(is(y, 'characteristics'))

  traverse.tree <- function(tree, level = NULL) {
    l <- lapply(names(tree),
                function(nodename) {
                  if ( is(tree[[nodename]], 'list') )
                    return(traverse.tree(tree[[nodename]],
                                         c(level, nodename)))

                  if ( verbose )
                    cat(sprintf('map: %s -> %s\n', paste(level, collapse = '.'),
                                                   nodename))

                  d <- x$dataparts(level, index = index)

                  if ( length(d) == 0 )
                    return(NA)

                  sapply(d, function(x) do.call(tree[[nodename]], unname(x)))
              })

    structure(l, names = names(tree))
  }

  structure(traverse.tree(y$map()),
            class = c('mapped.dataset', 'list'),
            name = attr(y, 'name'))
}



reduce <- function(x, y, ...) {
  UseMethod('reduce')
}


reduce.mapped.dataset <- function(x, y, verbose = TRUE, ...) {
  stopifnot(is(y, 'characteristics'))

  traverse.tree <- function(tree, level = NULL) {
    lapply(names(tree),
           function(nodename) {
             if ( is(tree[[nodename]], 'list') )
               return(traverse.tree(tree[[nodename]],
                                    c(level, nodename)))

             if ( verbose )
               cat(sprintf('reduce: %s\n', paste(c(level, nodename), collapse = '.')))

             f <- tree[[nodename]]

             if ( is.function(f) )
               x[[c(level, nodename)]] <<- f(x[[c(level, nodename)]])

             if ( is.null(f) )
               x[[c(level, nodename)]] <<- NULL

             if ( is(f, 'p') )
               x[[c(level, nodename)]] <<- do.call(f$fn,
                                                   lapply(f$args,
                                                          function(.) x[[.]]))
           })
  }

  traverse.tree(y$reduce())

  structure(x, class = c('reduced.dataset', class(x)),
            name = attr(y, 'name'))
}


