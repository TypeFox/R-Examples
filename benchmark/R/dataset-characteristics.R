#' @include proto.R
{}



DatasetCharacteristics <- proto::proto(expr = {
  name <- "Generic"

  requirements <- function(., ...) NULL

  map <- function(., ...) list()
  reduce <- function(., ...) list()


  pprint <- function(., ...) {
    cat(.$name, "characteristics\n")
  }

  psummary <- function(., ...) {
    chars <- .$characteristics(which = "reduce", flat = TRUE)

    .$print(...)
    cat(paste(" ", chars, collapse = "\n"), "\n")
  }

  characteristics <- function(., which = c("reduce", "map"), flat = TRUE, ...) {
    traverse.tree <- function(tree, level = NULL) {
      l <- lapply(names(tree),
                  function(nodename) {
                    if ( is.null(tree[[nodename]]) )
                      return(NULL)

                    if ( is(tree[[nodename]], "list") )
                      return(traverse.tree(tree[[nodename]],
                                           c(level, nodename)))

                    NA
                })

      structure(l, names = names(tree))
    }

    which <- match.arg(which)
    chars <- traverse.tree(do.call(which, list(), envir = .))

    if ( flat )
      names(unlist(chars))
    else
      chars
  }
})



### Definition helper functions: #####################################

o <- function(...) {
  fs <- list(...)
  function(...) Reduce(function(x, f) f(x), fs, ...)
}



p <- function(fn, args) {
  structure(list(fn = fn, args = args), class = 'p')
}



### Implementation -- StatLog characteristics: #######################


#' StatLog dataset characteristics
#'
#' Implementation of the StatLog project dataset characteristics.
#'
#' @usage
#'   StatlogCharacteristics
#'
#' @family dataset-characterization
#'
#' @references
#'   See \emph{Eugster et al. (2010)} in \code{citation("benchmark")}.
#'
#'   R. D. King, C. Feng and A. Sutherland. STATLOG: Comparison of
#'   classification algorithms on large real-world problems. Applied
#'   Artifical Intelligence, 9, 1995.
#'
#' @export
StatlogCharacteristics <- proto::proto(DatasetCharacteristics, expr = {
  name <- "Statlog"

  requirements <- function(., ...) {
    stopifnot(require(e1071))
    stopifnot(require(entropy))

    TRUE
  }

  map <- function(.) {
    map <- list()

    map$input <- list(n = nrow,
                      attr = ncol,
                      factor = list(attr = ncol,
                                    . = list(nlevels = nlevels,
                                             entropy = o(na.omit, as.integer,
                                                         entropy.empirical))),
                      numeric = list(attr = ncol,
                                     mac = mac,
                                     . = list(skewness = o(na.omit, skewness),
                                              kurtosis = o(na.omit, kurtosis))))

    map$response <- list(factor = list(. = list(cl = nlevels,
                                                entropy = o(na.omit, as.integer,
                                                            entropy.empirical))))

    map$input2response <- list(numeric2factor = list(fcc = fcc,
                                                     frac1 = frac1),
                               factor2factor = list(. = list(mi = mi)))
    map
  }

  reduce <- function(.) {
    reduce <- list()

    reduce$input <- list(n = identity,
                         attr = identity,
                         factor = list(attr = na0,
                                       . = list(bin = p(binary, list(c("input", "factor", ".", "nlevels"))),
                                                entropy = mean,
                                                nlevels = NULL)),
                         numeric = list(attr = na0,
                                        mac = mean,
                                        . = list(skewness = mean,
                                                 kurtosis = mean)))

    reduce$response <- list(factor = list(. = list(cl = identity,
                                                   entropy = identity)))

    reduce$input2response <- list(numeric2factor = list(fcc = identity,
                                                        frac1 = identity),
                                  factor2factor = list(. = list(mi = mean),
                                                       enattr = p(enattr, list(c("response", "factor", ".", "entropy"),
                                                                               c("input2response", "factor2factor", ".", "mi"))),
                                                       nsratio = p(nsratio, list(c("input", "factor", ".", "entropy"),
                                                                                 c("input2response", "factor2factor", ".", "mi")))))
    reduce
  }
})

StatlogCharacteristics <- structure(StatlogCharacteristics,
                                    class = c("characteristics",
                                              class(StatlogCharacteristics)))



### Implementation of needed characteristics: ########################

na0 <- function(x) {
  ifelse(is.na(x), 0, x)
}

enc <- function(x) {
  y <- matrix(0, nrow=length(x), ncol=nlevels(x))
  y[cbind(seq(length(x)), as.numeric(x))] <- 1
  y
}

mac <- function(x) {
  x <- as.matrix(x)

  if ( ncol(x) == 1 )
    return(NA)

  drop(sapply(seq(length = ncol(x)),
              function(i)
              sqrt(summary(lm(x[, i] ~ x[, -i]))$r.squared)))
}

fcc <- function(x, y) {
  x <- as.matrix(x)
  y <- unlist(y)
  max(cancor(x, enc(y))$cor)
}

frac1 <- function(x, y) {
  x <- as.matrix(x)
  y <- unlist(y)
  cor <- cancor(x, enc(y))$cor
  lambda <- cor^2 / (1 - cor^2)

  max(lambda) / sum(lambda)
}

mi <- function(x, y) {
  mi.plugin(cbind(as.integer(x),
                  as.integer(y)))
}

binary <- function(x, ...) {
  na0(sum(x == 2))
}

enattr <- function(x, y, ...) {
  x / y
}

nsratio <- function(x, y, ...) {
  (x - y) / y
}
