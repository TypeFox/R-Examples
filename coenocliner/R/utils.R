##' @title List parameters of species response models
##'
##' @description Returns the parameters of the indicated response model.
##'
##' @param model character; the species response model for which parameters will be listed
##'
##' @return A character vector of parameters. The species response model is returned as attribute \code{"model"}. Attribute \code{"onlyx"} is a logical vector indicating which, if any, of the parameters are intended to be supplied only once per species and not for both gradients.
##'
##' @author Gavin L. Simpson
##'
##' @export
##'
##' @seealso \code{\link{Gaussian}} and \code{\link{Beta}} for the species response model functions themselves.
##'
##' @examples
##' showParams("gaussian")
`showParams` <- function(model = c("gaussian", "beta")) {
    if (missing(model)) {
        stop("The species response model be provided.")
    }
    model <- match.arg(model)
    switch(model,
           gaussian = {
               params <- c("opt", "tol", "h")
               attr(params, "onlyx") <- c(rep(FALSE, 2), TRUE)
               attr(params, "model") <- "Gaussian"
           },
           beta = {
               params <- c("m", "r", "alpha", "gamma", "A0")
               attr(params, "onlyx") <- c(rep(FALSE, 4),TRUE)
               attr(params, "model") <- "Generalise Beta"
           })
    class(params) <- "modelParams"
    params
}

##' @export
`print.modelParams` <- function(x, ...) {
    onlyx <- attr(x, "onlyx")
    if (length(x)) {
        msg <- paste("Species response model:", attr(x, "model"))
        writeLines(strwrap(msg))
        writeLines(strwrap("Parameters:", prefix = "\n"))
        print(paste0(x, ifelse(onlyx, "*", "")), quote = FALSE, ...)
        if (any(onlyx)) {
            msg <- "Parameters marked with '*' are only supplied once"
            writeLines(strwrap(msg, prefix = "\n\t"))
        }
    }
    invisible(x)
}

##' @title Stack columns of a species coenocline simulation
##'
##' @description Stacks columns of a species coenocline simulation into long format suitable for use in statistical modeling or ggplot/lattice plots.
##'
##' @param x an object of class \code{"coenocline"}
##' @param ... arguments passed to other methods (not used).
##'
##' @importFrom utils stack
##'
##' @export
##'
##' @seealso \code{\link{coenocline}}
`stack.coenocline` <- function(x, ...) {
    X <- as.vector(x)
    vars <- colnames(x)
    if (is.null(vars)) {
        vars <- paste0("Spp", seq_len(NCOL(x)))
    }
    ind <- factor(rep(vars, each = NROW(x)))
    out <- cbind.data.frame(Species = ind, Abundance = X)
    out
}

##' @title Pretty Display of a Matrix or Data Frame
##'
##' @description Pretty display of the first \code{n} rows of a data frame or matrix-like object, with variables/columns that cannot fit on a single screen width removed.
##' @param x an R object for which \code{\link{head}} and \code{\link{as.data.frame}} methods exist.
##' @param n numeric; the number of rows to display
##' @param width numeric; the display width to assume when formatting the data frame. The default is given by \code{getOption("width")}
##' @param zapsmall logical; should values close to zero be zapped to zero? See \code{\link{zapsmall}} for details.
##' @return A \code{\link{format}}-ed data frame with \code{n} rows (or \code{n+1} rows if only a subset is shown) and as many columns/components that can be printed on a single screen width.
##' @author Gavin L. Simpson
##'
##' @rdname coenocliner-internal
##' @keywords internal
##' @importFrom utils head
##' @export
##'
##' @examples
##'
##' x <- seq(from = 4, to = 6, length = 30)
##' opt <- seq(4, 7, length = 100)
##' tol <- rep(0.25, 100)
##' h <- rep(20, 100)
##'
##' ## simulate
##' set.seed(1)
##' y <- coenocline(x, responseModel = "gaussian",
##'                 params = cbind(opt = opt, tol = tol, h = h),
##'                 countModel = "poisson")
##'
##' prettyHead(y)
##'
`prettyHead` <- function(x, n = 10, width = getOption("width"),
                         zapsmall = FALSE) {
    ## some meta data
    nr <- NROW(x)
    nc <- NCOL(x)
    ## truncate x to n rows using appropriate head() method
    df <- as.data.frame(head(x, n = n))
    if (zapsmall) {
        df <- zapsmall(df)
    }
    ## nullify rownames just incase there are any (shouldn't be...)
    rownames(df) <- NULL
    ## add some variable names; should make this user-configurable TODO
    colnames(df) <- paste0("Sp", seq_len(nc))

    ## format df converts everything to characters but stays in data
    ## frame format
    fdf <- format(df)
    ## this is a vector of character strings that represents
    ## what will get printed on a single line
    strings <- c(format(rownames(fdf))[[1]], unlist(fdf[1, ]))
    ## account for the rownames, hence initial empty string
    names <- c("", colnames(df))
    ## take length of elements of strings or names, whichever is largest
    ## these are the widths required to display each column of fdf
    widths <- pmax(nchar(strings), nchar(names))
    ## cumulative sum is total width up to & including the nth column
    cum <- cumsum(widths + 1)
    ## and we find which of these exceeds the current display width
    overfl <- cum[-1] > width           # ignore rownames

    if (all(overfl)) {                  # sanity check, always show 1 column
        overfl[1] <- FALSE
    }
    ## should need the above unless /very/ small width with huge counts

    ## loose columns we can't show in a single screen width
    trimDF <- fdf[, !overfl, drop = FALSE]

    ## dplyr adds continuation ellipsis to each column
    ## as counts tend to have narrower columns than dplyrs use-case
    ## add ":" as continuation chars in spirit of TeX \vdots
    ellipsis <- nr > n
    if (ellipsis) {
        dots <- rep(":", ncol(trimDF))
        trimDF <- rbind(trimDF, ":" = dots)
    }

    ## record number of species omitted
    omitted <- if (any(overfl)) {
        sum(overfl)
    } else {
        0
    }

    ## return object
    structure(list(df = trimDF, omitted = omitted, zapped = zapsmall),
              class = "prettyHead")
}

##' @export
`print.prettyHead` <- function(x, ...) {
    cat("\n")
    print(x$df, ...)
    cat("\n")
    if (x$omitted > 0) {
        writeLines(strwrap(paste("Counts for", x$omitted, "species not shown.")))
    }
    if (x$zapped) {
        writeLines(strwrap("(Values very close to zero were zapped)"))
    }
    cat("\n")
    invisible(x)
}

##' @title Return the First of Last Part of an Object
##'
##' @description Returns the first of last \code{n} rows of the simulated counts.
##'
##' @param x an object of class \code{"coenocline"}, the result of a call to \code{\link{coenocline}}.
##' @param ... additional arguments passed to other methods.
##'
##' @return A matrix consisting of the first or last \code{n} rows of the matrix of simulated counts.
##'
##' @author Gavin L. Simpson
##'
