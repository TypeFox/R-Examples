##' Create sequences for plotting
##'
##' \code{plotSeq} is a convenience for the creation of sequence
##' across the range of a variable.
##' By default, the length of the plotting
##' sequence will be equal to the length of the original sequence.
##' In that case, the only effect is to create an evenly-spaced
##' set of values. If \code{length.out} is specified, the user
##' determines the number of elements in plotSeq.
##'
##' The primary intended usage is for the creation of
##' plotting sequences of numeric variables.  It takes
##' a variable's range and the fills in evenly spaced steps.
##' If x is a factor variable, the levels will be returned.
##' Uses of this functionality are planned in the future.
##' @usage plotSeq(x, length.out = length(x))
##' @param x an R vector variable
##' @param length.out the number of elements in the desired plotting sequence.
##' @export plotSeq
##' @seealso \code{pretty}
##' @examples
##' #Create a quadratic regression
##'
##' stde <- 14
##' x <- rnorm(100, m = 50, s = 10)
##' y <- 0.2 - 02*x + 0.2*x^2 + stde*rnorm(100)
##' mod1 <- lm (y ~ poly(x, 2))
##'
##' plot(x, y, main="The Quadratic Regression")
##' seqx <- plotSeq(x, length.out = 10)
##' seqy <- predict(mod1, newdata = data.frame(x = seqx))
##' lines(seqx, seqy, col = "red")
##'
##' # Notice the bad result when a plotting sequence is
##' # not used.
##' plot(x, y, main = "Bad Plot Result")
##' seqy <- predict(mod1)
##' lines(x, seqy, col = "green")
##'
plotSeq <-
    function (x, length.out = length(x))
{
  if (is.numeric(x)) {
    xr <- range(x, na.rm = TRUE)
    pseq <- seq(xr[1], xr[2], length.out = length.out)
    return(pseq)
  } else {
    if (is.factor(x)) {
      pseq <- levels(x)
      return(pseq)
    } else {
      stop("plotSeq can only create plotting sequences for numeric or factor variables")
    }
  }
}
NULL


##' Select most frequently occurring values from numeric or categorical variables.
##'
##' The "n" most frequently occurring values are returned, sorted by
##' frequency of occurrence (in descending order). The names attribute
##' includes information about the percentage of cases that have the
##' indicated values.
##'
##' This is used by plotSlopes, plotCurves, and other "newdata" making
##' functions.
##' @param x A numeric or character variable
##' @param n The maximum number of values that may be returned.
##' @param pct Default = TRUE. Include percentage of responses within each category
##' @return A named vector.
##' @export
##' @author Paul E. Johnson <pauljohn@@ku.edu>
cutByTable <-
    function(x, n = 5, pct = TRUE)
{
    table1 <- table(x)
    n <- min(n, length(names(table1)))
    table1 <- sort(table1, decreasing = T)
    table1cut <- table1[1:n]
    tabNames <- names(table1cut)
    res <- if(is.factor(x)){
        as.factor(tabNames)
    } else if (is.numeric(x)){
        as.numeric(tabNames)
    } else {
        tabNames
    }
    if (pct) {
        freq <- 100*round(table1cut/sum(table1), 2)
        names(res) <- paste(tabNames," (",freq,"%)", sep="")
    } else {
        names(res) <- tabNames
    }
    if (is.factor(x)) res <- factor(res, levels = levels(x))
    res
}
NULL


##' Calculates the "center" quantiles, always including the median, when n is odd.
##'
##' If the numeric variable has fewer than 6 unique observed values,
##' this will send the data to \code{cutByTable}.  The default return
##' will find dividing points at three quantiles: c(0.25, 0.50, 0.75)
##' If n=4, the dividing points will be c(0.20, 0.40, 0.60, 0.80) If
##' n=5, c(0.0, 0.25, 0.50, 0.75, 1.0) Larger n that are odd will
##' include 0.5 and evenly spaced points out to proportions 0 and
##' 1.0. Larger n that is even will return evenly spaced points
##' calculated by R's \code{pretty} function.
##' @param x A numeric vector.
##' @param n The number of quantile points. See details.
##' @return A vector
##' @export
##' @importFrom stats quantile
##' @author Paul E. Johnson <pauljohn@@ku.edu>
cutByQuantile <-
    function(x, n = 3)
{
    uniqueVals <- unique(x)
    n <- as.integer(n) ##defend self against fractions
    if (length(uniqueVals) < 6) {
        qs <- cutByTable(x, n)
        invisible(qs)
    } else {
        cutVector <- if(n == 1) {
            c(0.50)
        } else if(n == 2) {
            c(0.25, 0.75)
        } else if(n == 3) {
            c(0.25, 0.50, 0.75)
        } else if(n == 4) {
            c(0.20, 0.40, 0.60, 0.80)
        } else if(n == 5) {
            c(0.00, 0.25, 0.50, 0.75, 1.0)
        } else if (n > 5) {
            if (n %% 2 == 0) {
                g <- 0.5 / n %/% 2
                c(seq(0, 0.5-g, by = g), 0.5, seq(0.5+g, 1.0, by = g))
            } else {
                pretty(c(0, 1), n = n-1)
            }
        }
        qs <- quantile(x, probs = cutVector, na.rm = TRUE)
        invisible(qs)
    }
}


##' Returns center values of x, the mean, mean-std.dev, mean+std.dev
##'
##' If the numeric variable has fewer than 6 unique observed values,
##' this will send the data to cutByTable.
##' @return A named vector
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @param x A numeric variable
##' @param n Should be an odd number 1, 3, 5, or 7. If 2 < n < 5,
##' values that divide the data at c(m-sd, m, m+sd) are returned. If n
##' > 4, the returned values are c(m-2sd, m-sd, m, m+sd, m+2sd).
##' @export
##' @importFrom stats sd
##' @examples
##' x <- rnorm(100, m = 100, s = 20)
##' cutBySD (x, n = 3)
##' cutBySD (x, n = 5)
cutBySD <-
    function(x, n = 3)
{
    uniqueVals <- unique(x)
    if (length(uniqueVals) < 6) {
        qs <- cutByTable(x, n)
        invisible(qs)
    } else {
        mx <- round(mean(x, na.rm = TRUE),2)
        sdx <- round(sd(x, na.rm = TRUE),2)
        if (n == 1) {
            qs <- c( mx )
            suffix <- c("(m)")
        } else if (n == 2) {
            qs <- c(mx - sdx, mx + sdx)
            suffix <- c("(m-sd)", "(m+sd)")
        } else if (n <= 4) {
            qs <- c(mx - sdx, mx, mx + sdx)
            suffix <- c("(m-sd)", "(m)", "(m+sd)")
        } else if (n <= 6){
            qs <- c(mx - 2*sdx, mx - sdx, mx, mx + sdx, mx + 2*sdx)
            suffix <- c("(m-2sd)","(m-sd)","(m)","(m+sd)","(m+2sd)")
        }  else {
            qs <- c(mx - 3*sdx, mx - 2*sdx, mx - sdx, mx, mx + sdx, mx + 2*sdx, mx + 3*sdx)
            suffix <- c("(m-3sd)", "(m-2sd)", "(m-sd)", "(m)", "(m+sd)", "(m+2sd)", "(m+3sd)")
        }
        names(qs) <-  paste(suffix)
        invisible(qs)
    }
}





##' Select focal values from an observed variable.
##'
##' This is a generic function with 2 methods, getFocal.default
##' handles numeric variables, while getFocal.factor handles
##' factors. No other methods have been planned for preparation.
##'
##' This is used in functions like \code{plotSlopes} or
##' \code{plotCurves}.
##'
##' @param x Required. A variable
##' @param ... Other arguments that will be passed to the
##' user-specified xvals function.
##' @return A vector.
##' @rdname getFocal
##' @export
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##'
getFocal <-
    function(x, ...)
{
    UseMethod("getFocal")
}

NULL

##' Select focal values from a numeric variable
##'
##' Many plotting functions need to select "focal" values from
##' a variable. There is a family of functions that are used
##' to do that.  User requests can be accepted in a number
##' of ways. Numeric and character variables will be treated
##' differently. Please see details.
##'
##' If \code{xvals} is not provided, a default divider for numeric
##' variables will be the algorithm "quantile". The divider algorithms
##' provided with rockchalk are c("quantile", "std.dev.", "table",
##' "seq").  \code{xvals} can also be the name of a user-supplied
##' function, such as R's \code{pretty()}. If the user supplies a
##' vector of possible values, that selection will be checked to make
##' sure all elements are within a slightly expanded range of
##' \code{x}. If a value out of range is requested, a warning is
##' issued. Maybe that should be an outright error?
##' 
##' With factor variables, \code{xvals} is generally not used because
##' the only implemented divider algorithm is "table" (see
##' \code{cutByTable}), which selects the \code{n} most frequently
##' observed values. That is the default algorithm. It is legal to
##' specify xvals = "table", but there is no point in doing
##' that. However, xvals may take two other formats. It may be a
##' user-specified function that can select levels values from
##' \code{x} or it may be a vector of labels (or, names of
##' levels). The purpose of the latter is to check that the requested
##' levels are actually present in the supplied data vector
##' \code{x}. If the levels specified are not in the observed
##' variable, then this function stops with an error message.
##'
##' @param xvals A function name (either "quantile", "std.dev.",
##' "table", or "seq") or a user-supplied function that can receive
##' x and return a selection of values.
##' @param n Number of values to be selected.
##' @param pct Default TRUE. Include percentage of observed cases
##' in variable name? (used in legends)
##' @return A named vector of values.
##' @export
##' @method getFocal default
##' @rdname getFocal
##' @examples
##' x <- rnorm(100)
##' getFocal(x)
##' getFocal(x, xvals = "quantile")
##' getFocal(x, xvals = "quantile", n = 5)
##' getFocal(x, xvals = "std.dev")
##' getFocal(x, xvals = "std.dev", n = 5)
##' getFocal(x, xvals = c(-1000, 0.2, 0,5))
##'
getFocal.default <- function(x, xvals = NULL, n = 3, pct = TRUE, ...)
{
    xRange <- magRange(range(x, na.rm = TRUE), 1.1)

    if (is.null(xvals))
        return(xfocal <- cutByQuantile(x, n))

    if (is.numeric(xvals)) {
        if ((xvals > max(xRange)) || (xvals < min(xRange))){
            warning("values requested out of observed range in getFocal")
        }
        xvals <- sort(xvals)
        return(xvals)
    }

    if (is.character(xvals)) {
        xvals <- match.arg(tolower(xvals),
                           c("quantile", "std.dev.","table", "seq"))
        xfocal <- switch(xvals,
                         table = cutByTable(x, n, pct = pct),
                         quantile = cutByQuantile(x, n),
                         "std.dev." = cutBySD(x, n),
                         seq = plotSeq(x, n),
                         stop("unknown 'xvals' algorithm in getFocal"))
        return(xfocal)
    }

    if (is.function(xvals)) {
        return(xvals(x, n, ...))
    }

    stop("getFocal received unexpected input for xvals")
}
NULL


##' @export
##' @method getFocal factor
##' @rdname getFocal
##' @examples
##' x <- factor(c("A","B","A","B","C","D","D","D"))
##' getFocal(x)
##' getFocal(x, n = 2)
##'
getFocal.factor <- function(x, xvals = NULL, n = 3, pct = TRUE, ...)
{
    if (is.null(xvals)) {
        xvals <- cutByTable(x, n, pct = pct)
        return(xvals)
    }
    if (is.vector(xvals)) {
       x <- factor(x) ## drop unused levels
       if (!all(xvals %in% levels(x))) stop("xvals includes non-observed levels of x")
       return(xvals)
    }
    if (is.function(xvals)) {
        xvals <- xvals(x, n, ...)
        return(xvals)
    }

    if (is.character(x)) {
        xvals <- match.arg(tolower(xvals),
                           c("table"))
        xvals <- switch(xvals,
                        table = cutByTable(x, n),
                        stop("Sorry, only known algorithm for factors is 'table'"))
        return(xvals)
    }

    if (is.function(xvals)) {
        return(xvals(x, n, ...))
    }

    stop("getFocal received unexpected input xvals")
}


##' @export
##' @method getFocal character
##' @rdname getFocal
##' @examples
##' x <- c("A","B","A","B","C","D","D","D")
##' getFocal(x)
##' getFocal(x, n = 2)
##'
getFocal.character <- function(x, xvals = NULL, n = 3, pct = TRUE, ...)
{
    if (is.null(xvals)) {
        xvals <- cutByTable(x, n, pct = pct)
        return(xvals)
    }
    if (is.vector(xvals)) {
       if (!all(xvals %in% unique(x))) stop("xvals includes non-observed values of x")
       return(xvals)
    }
    if (is.function(xvals)) {
        xvals <- xvals(x, n, ...)
        return(xvals)
    }

    xvals <- match.arg(tolower(xvals),
                       c("table"))
    xvals <- switch(xvals,
                    table = cutByTable(x, n, pct = pct),
                    stop("Sorry, only known algorithm for character variables is 'table'"))
    xvals
}



##' Pad with 0's.
##'
##' Sometimes we receive this c(1, 22, 131) and we need character
##' variables of the same size, such as c("001", "022", "131").  This
##' happens if a user has mistakenly converted a zip code (US regional
##' identifier) like "00231" to a number. This function converts the
##' number back to a 0 padded string.
##'
##' This works differently if the number provided is an integer, or a
##' character string.  Integers are left padded with the character
##' "0".  A character string will be left-padded with blanks.
##' @param x a numeric  variable.
##' @return A character string vector padded with 0's
##' @export
##' @author Paul Johnson <pauljohn@@ku.edu>
##' @examples
##' x <- c(1 , 11, 22, 121, 14141, 31)
##' (xpad <- padW0(x))
##' x <- rpois(7, lambda = 11)
##' (xpad <- padW0(x))
##' x <- c("Alabama", "Iowa", "Washington")
padW0 <- function(x){
    vtype <- "d"
    if (is.character(x)) vtype <- "s"
    sprintf(paste("%0", max(nchar(as.character(x))), vtype, sep=""), x)
}


##' A way of checking if a string is a valid file name.
##'
##' A copy from R's grDevices:::checkIntFormat because it is not exported there
##' @param s An integer
##' @return logical: TRUE or FALSE
##' @author R Core Development Team
checkIntFormat <- function(s) {
    s <- gsub("%%", "", s)
    if (length(grep("%", s)) == 0L) 
        return(TRUE)
    s <- sub("%[#0 ,+-]*[0-9.]*[diouxX]", "", s)
    length(grep("%", s)) == 0L
}


##' Create a uniquely named directory. Appends number & optionally date to directory
##' name.
##' 
##' Checks if the requested directory exists. If so, will create new directory
##' name. My favorite method is to have the target directory with a date-based
##' subdirectory, but set usedate as FALSE if you don't like that. Arguments
##' showWarnings, recursive, and mode are passed along to R's dir.create, which
##' does the actual work here.
##'
##' Default response to dir = "../output/" fixes the directory name like this,
##' "../output/20151118-1/" because usedate is assumed TRUE. If usedate = FALSE,
##' then output names will be like "../output-1/", "../output-2/", and so forth.
##' @param path A character string for the base name of the directory.
##' @param usedate TRUE or FALSE: Insert YYYYMMDD information?
##' @param showWarnings default TRUE. Show warnings? Will be passed on to dir.create
##' @param recursive default TRUE. Will be passed on to dir.create
##' @param mode  Default permissions on unix-alike systems. Will be passed on to dir.create
##' @export
##' @return a character string with the directory name
##' @author Paul E Johnson <pauljohn@@ku.edu>
dir.create.unique <- function(path, usedate = TRUE,
                              showWarnings = TRUE,
                              recursive = TRUE,
                              mode = "0777"){
    dts <- function(name) gsub("/$", "", name)
    if (!file.exists(dts(path))){
        dir.create(path, recursive = recursive, showWarnings = showWarnings, mode = mode)
        return(path)
    } else {
        today <- format(Sys.time(), "%Y%m%d")
        j <- 1
        ocandidate <- paste0(dts(path), ifelse(usedate, paste0("/", today), ""), "-",  j, "/")
        while (file.exists(dts(ocandidate))) {
            j <- j + 1
            ocandidate <- paste0(dts(path), ifelse(usedate, paste0("/", today), ""), "-",  j, "/")
        }
        path <- ocandidate
        dir.create(path, recursive = recursive, showWarnings = showWarnings, mode = mode)
        return(path)
    }
}

