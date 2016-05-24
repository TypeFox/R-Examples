##' @title An \code{expand.grid}-like function that repeats sets of
##' vectors for every value in a reference vector.
##'
##' @description The values of \code{x} are repeated for each combination
##' of elements in the vectors supplied via \code{...}, with the first
##' elements of each vector in \code{...} being taken as a set, the
##' second elements as another set, and so on. \code{x} is repeated for
##' each of these sets.
##'
##' @param x numeric; vector of data points which are to be replicated
##' for each of the sets of vectors supplied to \code{...}.
##' @param ... additional vector arguments to be expanded to the correct
##' length. These are taken to be a set of values to be replicated for
##' each of the elements of \code{x}.
##'
##' @return a matrix of replicated vectors, with column names for \code{x}
##' and named arguments passed as \code{...}.
##'
##' @author Gavin L. Simpson
##'
##' @export
##'
##' @keywords utilities
##'
##' @references Minchin P.R. (1987) Simulation of multidimensional
##' community patterns: towards a comprehensive model. \emph{Vegetatio}
##' \strong{71}, 145--156.
##'
##' @examples
##' # Recreate Fig. 2 of Minchin (1987)
##' # Parameters for each of 6 six species
##' A0 <- c(5,4,7,5,9,8) * 10
##' m <- c(25,85,10,60,45,60)
##' r <- c(3,3,4,4,6,5) * 10
##' alpha <- c(0.1,1,2,4,1.5,1)
##' gamma <- c(0.1,1,2,4,0.5,4)
##' # Gradient locations
##' x <- 1:100
##'
##' # expand parameter set
##' pars <- expand(x, m = m, A0 = A0, r = r, alpha = alpha,
##'                gamma = gamma)
##' head(pars)
`expand` <- function(x, ...) {
    ## if `...` is length 1 & matrix use expandMatrix()
    ## otherwise, expandList
    dots <- list(...)
    out <- if (length(dots) == 1L && is.matrix(dots[[1]])) {
        expandMatrix(x, params = dots[[1]])
    } else {
        args <- vector("list", length = length(dots) + 1)
        args[[1]] <- x
        args[-1] <- dots
        names(args) <- c("x", names(dots))
        do.call("expandList", args)
    }
    out
}

`expandList` <- function(x, ...) {
    expandFun <- function(x, r1, r2) {
        x[rep.int(rep.int(seq_len(r2), rep.int(r1, r2)), 1L)]
    }
    dots <- list(...)
    nams <- names(dots)
    ## nams could be NULL or have zero length values
    namsD <- paste0("Var", seq_along(dots))
    if (is.null(nams)) {
        nams <- namsD
    } else if (any(want <- nzchar(nams))){
        namsD[want] <- nams[want]
    }
    nx <- length(x)
    n1 <- length(dots[[1]])
    orep <- nx * n1
    x <- x[rep.int(rep.int(seq_len(nx), rep.int(1L, nx)), orep / nx)]
    other <- vapply(dots, FUN = expandFun,
                    FUN.VALUE = numeric(length = length(x)),
                    r1 = nx, r2 = n1)
    out <- cbind(x, other)
    colnames(out)[-1] <- namsD
    out
}

`expandMatrix` <- function(x, params) {
    nx <- length(x)
    np <- NROW(params)
    orep <- nx * np
    x <- x[rep.int(rep.int(seq_len(nx), rep.int(1L, nx)), orep / nx)]
    other <- params[rep.int(rep.int(seq_len(np), rep.int(nx, np)), 1L), ]
    out <- cbind(x, other)
    colnames(out) <- c("x", colnames(params))
    out
}
