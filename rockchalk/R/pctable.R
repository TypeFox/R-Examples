##' Creates a cross tabulation with counts and percentages
##'
##' @description
##' This function is pronounced "presentable"!  The original purpose
##' was to create a particular kind of cross tabulation that I ask for
##' in class: counts with column percentages. Requests from users have
##' caused a bit more generality to be built into the function. Now,
##' optionally, it will provide row percents. This is a generic function.
##' Most users will find the formula method most convenient. Use the
##' colpct and rowpct arguments to indicate if column or row percentages
##' are desired.
##'
##' I suggest most users will use the formula method for this. Running
##' a command like this will, generally, do the right thing:
##'
##' \code{tab <- pctable(y ~ x, data = dat)}
##'
##' There is also a method that will work with characters representing
##' variable names.
##' 
##' \code{tab <- pctable("y", "x", data = dat)}
##'  
##' Running the function should write a table in the output console,
##' but it also creates an object (\code{tab}). That object
##' can be displayed in a number of ways.
##'
##' A summary method is provided, so one could look at different
##' representations of the same table.
##'
##' \code{summary(tab, rowpct = TRUE, colpct = FALSE)}
##'
##' or
##'
##' \code{summary(tab, rowpct = TRUE, colpct = TRUE)}
##'
##' Tables that include only row or column percentages will be
##' compatible with the html and latex exporters in the
##' excellent \code{tables} package.
##' 
##' @details
##' Please bear in mind the following. The output object is a
##' list of tables of partial information, which are then assembled in
##' various ways by the print method (print.pctable). A lovely table
##' will appear on the screen, but the thing tab itself has more
##' information and a less beautiful structure. In order to see the
##' individual pieces, run tab[[1]], or tab[[2]], or tab[[3]].
##'
##' The object tab is of class \code{pctable}. A print method is supplied.
##' For any pctable object, it is possible to run follow-ups like
##'
##' print(tab, rowpct = TRUE, colpct = FALSE)
##' 
##' The method \code{print.pctable(tab)} assembles the object into (my
##' opinion of) a presentable form. The print method has argumnets
##' \code{rowpct} and \code{colpct} that determine which percentages
##' are included in the presentation.
##' 
##' @name pctable
NULL

##' @param rv A row variable name
##' @param ... Other arguments. So far, the most likely additional
##' arguments are to be passed along to the table function, such as
##' "exclude", "useNA", or "dnn" (which will override the rvlab and
##' cvlab arguments provided by some methods). Some methods will
##' also pass along these arguments to model.frame, "subset" 
##' "xlev", "na.action", "drop.unused.levels".
##' @seealso \code{\link[tables]{tabular}} and the CrossTable function
##' in \code{gmodels} package.
##' @rdname pctable
##' @export
##' @return A list with tables (count, column percent, row percent) as
##' well as a copy of the call.
##' @author Paul Johnson <pauljohn@@ku.edu>
pctable <- function(rv, ...)
{
    UseMethod("pctable")
}
NULL




##' The method pctable.default is the calculator, I don't expect
##' many users will need to call it directly.
##' 
##' @param cv Column variable
##' @param rvlab Optional: row variable label
##' @param cvlab Optional: col variable label
##' @param colpct Default TRUE: are column percentags desired in the
##' presentation of this result?
##' @param rowpct Default FALSE: are row percentages desired in the
##' presentation of this result
##' @param rounded Default FALSE, rounds to 10's for privacy purposes.
##' @rdname pctable
##' @method pctable default
##' @rdname pctable
##' @export
pctable.default <- function(rv, cv,
                            rvlab = NULL, cvlab = NULL,
                            colpct = TRUE, rowpct = FALSE,
                            rounded = FALSE, ...)
{
    rvlabel <- if (!missing(rv)) deparse(substitute(rv))
    cvlabel <- if (!missing(cv)) deparse(substitute(cv))
    rvlab <- if (is.null(rvlab)) rvlabel else rvlab
    cvlab <- if (is.null(cvlab)) cvlabel else cvlab

    dots <- list(...)
    dotnames <- names(dots)
    ## altargs <- list()
   
    ## if ("dnn" %in% dotnames) altargs$dnn <- dots[["dnn"]]
    ## if ("deparse.level" %in% dotnames)
    ##     altargs$deparse.level <- dots[["deparse.level"]]
    
    tableargs <- list(rv, cv, dnn = c(rvlab, cvlab))
    newargs <- modifyList(tableargs, dots, keep.null = TRUE)

    t1 <- do.call("table", newargs)
    rownames(t1)[is.na(rownames(t1))] <- "NA" ## symbol to letters
    colnames(t1)[is.na(colnames(t1))] <- "NA"
    if (rounded) t1 <- round(t1, -1)
    t2 <- addmargins(t1, c(1,2))
    t1cpct <- round(100*prop.table(t1, 2), 1)
    t1rpct <- round(100*prop.table(t1, 1), 1)
    t1cpct <- apply(t1cpct, c(1,2), function(x) gsub("NaN", "", x))
    t1rpct <- apply(t1rpct, c(1,2), function(x) gsub("NaN", "", x))
    res <- list("count" = t2, "cpct" = t1cpct, "rpct" = t1rpct,
                call = match.call())
    class(res) <- "pctable"
    print.pctable(res, colpct = colpct, rowpct = rowpct)
    invisible(res)
}
NULL




##' Creates a cross tabulation with counts and column percents
##'
##' The formula method is the recommended method for users. Run
##' \code{pctable(myrow ~ mycol, data = dat)}. In an earlier version,
##' I gave different advice, so please adjust your usage.
##' 
##' @param formula A two sided formula.  
##' @param data A data frame.
##' @examples
##' dat <- data.frame(x = gl(4, 25),
##'                   y = sample(c("A", "B", "C", "D", "E"), 100, replace= TRUE))
##' pctable(y ~ x, dat)
##' pctable(y ~ x, dat, exclude = NULL)
##' pctable(y ~ x, dat, rvlab = "My Outcome Var", cvlab = "My Columns")
##' pctable(y ~ x, dat, rowpct = TRUE, colpct = FALSE)
##' pctable(y ~ x, dat, rowpct = TRUE, colpct = TRUE)
##' pctable(y ~ x, dat, rowpct = TRUE, colpct = TRUE, exclude = NULL)
##' tab <- pctable(y ~ x, dat, rvlab = "Outcome", cvlab = "Predictor")
##' @rdname pctable
##' @method pctable formula
##' @export
pctable.formula <- function(formula, data = NULL,  rvlab = NULL,
                            cvlab = NULL, colpct = TRUE, rowpct = FALSE,
                            rounded = FALSE, ...)
    
{
    if (missing(formula) || (length(formula) != 3L))
        stop("pctable requires a two sided formula")
    dots <- list(...)
    dotnames <- names(dots)
    
    mt <- terms(formula, data = data)
    if (attr(mt, "response") == 0L) stop("response variable is required")
    mf <- match.call(expand.dots = TRUE)
    mfnames <- c("formula", "data", "subset", "xlev", "na.action", "drop.unused.levels")
    keepers <- match(mfnames, names(mf), 0L)
    mf <- mf[c(1L, keepers)]
    
    ## mf$drop.unused.levels <- FALSE
    if (!"na.action" %in% dotnames) mf$na.action <- na.pass
    mf[[1L]] <- quote(stats::model.frame)

    ## remove used arguments from dots, otherwise errors happen
    ## when unexpected arguments pass through. Don't know why
    for (i in c("subset", "xlev", "na.action", "drop.unused.levels")) dots[[i]] <- NULL
        
    mf <- eval(mf, parent.frame())
    mfnames <- names(mf)
    response <- attr(attr(mf, "terms"), "response")
    ## response is column 1
    rvname <- mfnames[response]
    cvname <- mfnames[-response][1] ##just take 2?
    rvlab <- if (missing(rvlab)) rvname else rvlab
    cvlab <- if (missing(cvlab)) cvname else cvlab

    arglist <- list(rv = mf[[rvname]], cv = mf[[cvname]],
                    rvlab = rvlab, cvlab = cvlab,
                    colpct = colpct, rowpct = rowpct,
                    rounded = rounded)
    arglist <- modifyList(arglist, dots, keep.null = TRUE)
    ##keep.null needed because exclude = NULL can be a valid
    ## (meaningful) argument to table.
    
    res <- do.call(pctable.default, arglist)
    invisible(res)
}
NULL





##' Method for variable names as character strings 
##'
##' The character method exists only for variety.  It accepts
##' character strings rather than a formula to define the columns that
##' should be plotted.  The method used most often for most users should
##' be the formula method.
##' 
##' When using character arguments, the row variable rv rowvar must be
##' a quoted string if the user intends the method pctable.character
##' to be dispatched. The column variable cv may be a string or just a
##' variable name (which this method will coerce to a string).
##' @examples
##' dat <- data.frame(x1 = gl(4, 25, labels = c("Good", "Bad", "Ugly", "Indiff")),
##'                 x2 = gl(5, 20, labels = c("Denver", "Cincy", "Baltimore", "NY", "LA")), 
##'                 y = sample(c("A", "B", "C", "D", "E"), 100, replace= TRUE))
##' tab <- pctable(y ~ x1, data = dat, rvlab = "my row label",
##'     subset = dat$x1 %in% c("Good", "Bad"),
##'     drop.unused.levels = TRUE)
##' tab <- pctable(y ~ x1, data = dat, rvlab = "my row label",
##'     subset = dat$x1 %in% c("Good", "Bad"))
##' pctable("y", "x1", dat)
##' pctable("y", x1, dat)
##' tab <- pctable(y ~ x2, data = dat, rvlab = "A Giant Variable")
##' summary(tab, rowpct = TRUE, colpct = FALSE)
##' tabsum <- summary(tab)
##' \donttest{
##' ## if user has tables package, can push out to latex or html
##' if (require(tables) & require(Hmisc)){
##'     tabsumtab <- tables::as.tabular(tabsum)
##'     Hmisc::html(tabsumtab)
##'     fn <- tempfile(pattern = "file", tmpdir = tempdir(), 
##'             fileext = ".html")
##'     Hmisc::html(tabsumtab, file = fn)
##'     print(paste("The file saved was named", fn, "go get it."))
##'     if (interactive()) browseURL(fn)
##'     ## go get the fn file if you want to import it in document
##'     ## Now LaTeX output
##'     ## have to escape the percent signs
##'     tabsumtab <- apply(tabsumtab, 1:2, function(x) {gsub("%", "\\\\%", x) })
##'     fn2 <- tempfile(pattern = "file", tmpdir = tempdir(), 
##'                    fileext = ".tex")
##'     Hmisc::latex(tabsumtab, file = fn2)
##'     print(paste("The file saved was named", fn2, "go get it."))
##' }
##' }
##' @rdname pctable
##' @method pctable character
##' @export
pctable.character <- function(rv, cv, data = NULL, rvlab = NULL,
                              cvlab = NULL, colpct = TRUE,
                              rowpct = FALSE,
                              rounded = FALSE, ...)
{
    if (missing(data) || !is.data.frame(data)) stop("pctable requires a data frame")
    cv <- as.character(substitute(cv))[1L]
    
    rvlab <- if (missing(rvlab)) rv else rvlab
    cvlab <- if (missing(cvlab)) cv else cvlab
    res <- pctable.formula(formula(paste(rv, " ~ ", cv)), data = data,
                           rvlab = rvlab, cvlab = cvlab, colpct = colpct,
                           rowpct = rowpct, rounded = rounded, ...)
  
    invisible(res)
}
NULL


##' Extract presentation from a pctable object
##'
##' Creates a column and/or row percent display of a pctable
##' result
##' @param object A pctable object 
##' @param colpct Default TRUE: should column percents be included
##' @param rowpct Default FALSE: should row percents be included
##' @param ... Other arguments, currently unused
##' @return An object of class summary.pctable
##' @author Paul Johnson <pauljohn@@ku.edu>
##' @method summary pctable
##' @export
summary.pctable <- function(object, ..., colpct = TRUE, rowpct = FALSE)
{
    colpct <- if (missing(colpct)) object$call[["colpct"]] else colpct 
    rowpct <- if (missing(rowpct)) object$call[["rowpct"]] else rowpct
    
    count <- object[["count"]]
   
    t3 <- count
    attr(t3, which = "colpct") <- colpct
    attr(t3, which = "rowpct") <- rowpct
    class(t3) <- c("summary.pctable", "table")
    if (colpct && !rowpct) {
        cpct <- object[["cpct"]]
        for(j in rownames(cpct)){
            for(k in colnames(cpct)){
                if(j != "" & k != ""){
                    t3[j, k] <- paste0(count[j, k], "(", cpct[j, k], "%)")
                }
            }
        }
        return(t3)
    }
    
    ## rowpct == TRUE else would have returned
    rpct <- object[["rpct"]]
    for(j in rownames(rpct)){
        for(k in colnames(rpct)){
            if (j != "" & k != ""){
                t3[j, k] <- paste0(count[j, k], "(", rpct[j, k], "%)")
            }
        }
    }
    
    if (!colpct) {
        return(t3)
    } else {
        cpct <- object[["cpct"]]
        t4 <- array("", dim = c(1, 1) + c(2,1)*dim(object$cpct))
        t4[seq(1, NROW(t4), 2), ] <- t3
        rownames(t4)[seq(1, NROW(t4), 2)] <- rownames(t3)
        rownames(t4)[is.na(rownames(t4))] <- "" 
        colnames(t4) <- colnames(t3)
        for(j in rownames(object[["cpct"]])) {
            for(k in colnames(object[["cpct"]])){
                if(j != "" & k != "") {
                    t4[1 + which(rownames(t4) == j) ,k] <- paste0(cpct[j, k], "%")
                }
            }
        }
        t4 <- as.table(t4)
        names(dimnames(t4)) <- names(dimnames(count))
        attr(t4, which = "colpct") <- colpct
        attr(t4, which = "rowpct") <- rowpct
        class(t4) <- c("summary.pctable", "table")
        return(t4)
    }
}

##' print method for summary.pctable objects 
##'
##' prints pctab objects. Needed only to deal properly with quotes
##' 
##' @param x a summary.pctable object
##' @param ... Other arguments to print method
##' @return Nothing is returned
##' @author Paul Johnson <pauljohn@@ku.edu>
##' @method print summary.pctable
##' @export
print.summary.pctable <- function(x, ...){
    colpct <- attr(x, "colpct")
    rowpct <- attr(x, "rowpct")
    
    if (colpct && !rowpct) {
        cat("Count (column %)\n")
    } else if (!colpct && rowpct) {
        cat("Count (row %)\n")
    } else {
        cat("Count (row %)\n")
        cat("column %\n")
    }
    NextMethod(generic = "print", object = x, quote = FALSE, ...)
    ## Error in NextMethod(generic = "print", object = x, quote = FALSE, ...) :
    ## no calling generic was found: was a method called directly?
    ## Calls: pctable ... print.pctable -> print.summary.pctable -> NextMethod
    ## Execution halted
    ## print.table(x, quote = FALSE, ...)
}


##' Display pctable objects
##'
##' This is not very fancy. Note that the saved pctable object
##' has the information inside it that is required to write both
##' column and row percentages. The arguments colpct and rowpct
##' are used to ask for the two types.
##' 
##' @param x A pctable object
##' @param colpct Default TRUE: include column percentages?
##' @param rowpct Default FALSE: include row percentages?
##' @param ... Other arguments passed through to print method
##' @return A table object for the final printed table.
##' @author Paul Johnson <pauljohn@@ku.edu>
##' @method print pctable
##' @export
print.pctable <- function(x, colpct = TRUE, rowpct = FALSE, ...)
{
    colpct <- if (missing(colpct)) x$call[["colpct"]] else colpct 
    rowpct <- if (missing(rowpct)) x$call[["rowpct"]] else rowpct
    tab <- summary(x, colpct = colpct, rowpct = rowpct)
    print(tab, ...)
    invisible(tab)
}
NULL
