#' @name print
#' @title \code{print} methods
#' 
#' @include autoplotTen.R
#' 
#' @details Prints a \code{ten} object with 'nice' formatting.
#'  \cr
#' Options may be set for a session using e.g.
#'  \cr
#' options(survMisc.nColSP=4L)
#'  \cr
#' It is similar to the behavior of \code{print.data.table} but
#' has additional arguments controlling the number of columns
#' sent to the terminal.
#'
#' @param x An object of class \code{ten}.
#' @param ... Additional arguments (not implemented).
#'  \cr \cr
#' \bold{--print.ten}
#' @param maxRow Maximum number of rows to print.
#'  \cr
#' If \code{nrow(x) > maxRow}, just the first and last
#' \code{nRowP} (below) are printed.
#'  \cr
#' The default value is that used by \code{data.table}.
#' @param nRowP \bold{N}umber of rows to \bold{p}rint from
#'  the start and end of the object. Used if \code{nrow(x) > maxRow}.
#' @param pRowNames Print row names? 
#'  \cr
#' Default is \code{TRUE}.
#' @param maxCol Maximum number of columns to print.
#'  \cr
#' If \code{ncol(x) > maxCol}, just the first \code{nColSP}
#' and last \code{maxCol - nColSP} columns are printed.
#' @param nColSP \bold{N}umber of \bold{col}umns to \bold{p}rint from
#'  the \bold{s}tart of the object. Used if Used if \code{ncol(x) > maxCol}.
#' @param sigDig \bold{Sig}nificant \bold{dig}its. This is passed as an argument to
#'  \cr
#' ?signif
#'  \cr
#' when preparing the object for printing.
#'  \cr \cr
#' \bold{--print.tableAndPlot} and \bold{print.tableAndPlot}
#' @param hideTabLeg Hide table legend.
#' @param tabHeight Table height (relative to whole plot).
#'  \cr \cr
#' \bold{--print.COV}
#' @param n Similar to \code{n} from e.g. 
#'  \cr
#' ?utils::head
#'  \cr \cr
#' \bold{--print.lrt}
#' @param dist Which distribution to use for the statistics
#' when printing.
#'  \cr
#' Default (\code{dist="n"}) prints \eqn{Z} and \eqn{p} values
#' based on the normal distribution. 
#'  \cr
#' If \code{dist="c"}, gives values based on the 
#' \eqn{\chi^2}{chi-squared} distribution. 
#'  \cr
#' The results are the same. The default value is typically 
#' easier to read. Both options are given for completeness.
#' 
#' @return A printed representation of the object
#' is send to the terminal as a \emph{side effect} of
#' calling the function.
#'  \cr
#' The return value cannot be \code{assign}ed.
#' 
#' @author Chris Dardis. Based on existing work by Brian Diggs.
#' 
#' @seealso For \code{print.ten}:
#' @seealso data.table:::print.data.table
#' @seealso ?stats::printCoefmat
#' @seealso options()$datatable.print.nrows
#' @seealso sapply(c("datatable.print.nrows", "datatable.print.topn"), getOption)
#' 
#' @note
#' All numeric arguments to the function must be supplied as integers.
#'
#' @rdname print
#' @method print ten
#' @aliases print.ten
#' @export
#'
#' @examples
#' set.seed(1)
#' (x <- data.table::data.table(matrix(rnorm(1800), ncol=15, nrow=120)))
#' data.table::setattr(x, "class", c("ten", class(x)))
#' p1 <- print(x)
#' stopifnot(is.null(p1))
#' x[1:80, ]
#' x[0, ]
#' (data.table::set(x, j=seq.int(ncol(x)), value=NULL))
#' 
print.ten <- function(x, ...,
                      maxRow=getOption("datatable.print.nrows", 50L),
                      nRowP=getOption("datatable.print.topn", 5L),
                      pRowNames=TRUE,
                      maxCol=getOption("survMisc.maxCol", 8L),
                      nColSP=getOption("survMisc.nColSP", 7L),
                      sigDig=getOption("survMisc.sigDig", 2L)){
    if (nrow(x)==0L) {
        if (length(x)==0L) {
            cat("Null 'ten' object (0 rows and 0 cols)\n")
        } else {
            cat("Empty 'ten' object (0 rows) of ", length(x),
                " col", if (length(x) > 1L)
                "s", ": ", paste(head(names(x), 6), collapse = ", "),
                if (ncol(x) > 6)
                "...", "\n", sep = "")
        }
        return(invisible())
    }
    stopifnot(all(sapply(list(maxRow, nRowP, maxCol, nColSP),
                         as.integer)))
    stopifnot(maxRow > nRowP)
    ## lCol1 = last columns; needs to be at least one
    stopifnot((lCol1 <- maxCol - nColSP) >= 1)
    if (nrow(x) > maxRow) {
        toPrint1 <- rbind(head(x, nRowP),
                          tail(x, nRowP))
        ## row names
        rn1 <- c(seq_len(nRowP),
                 seq.int(to=nrow(x), length.out=nRowP))
        rowDots1 <- TRUE
    } else {
        toPrint1 <- x
        rn1 <- seq_len(nrow(x))
        rowDots1 <- FALSE
    }
    if (ncol(x) > (nColSP + lCol1 + 1L)) {
        toPrint1 <- cbind(
            toPrint1[, seq.int(nColSP), with=FALSE],
            toPrint1[, seq.int(to=ncol(x), length.out=lCol1), with=FALSE])
        colDots1 <- TRUE
    } else {
        colDots1 <- FALSE
    }
    toPrint1 <- do.call("cbind",
                        lapply(toPrint1,
                               function(col) signif(col,
                                                    digits=sigDig)))
    if (pRowNames) {
        rownames(toPrint1) <- paste(format(rn1, right = TRUE),
                                    ":",
                                    sep = "")
    } else {
        rownames(toPrint1) <- rep.int("", nrow(x))
    }
    if (rowDots1) {
        toPrint1 <- rbind(head(toPrint1, nRowP),
                          "---" = "",
                           tail(toPrint1, nRowP))
        rownames(toPrint1) <- format(rownames(toPrint1),
                                     justify="right")
    }
    if (colDots1) {
        toPrint1 <- cbind(
            toPrint1[, seq.int(nColSP), drop=FALSE],
            rep("", nrow(toPrint1)),
            toPrint1[, seq.int(to=ncol(toPrint1), length.out=lCol1),
                     drop=FALSE])
        colnames(toPrint1)[colnames(toPrint1)==""] <- " ---"
    }
    if (!rowDots1) {
        toPrint1 <- rbind(toPrint1,
                          matrix(colnames(toPrint1), nrow=1L))
    }
    print(toPrint1, right=TRUE, quote=FALSE)
    return(invisible())
}
#'
#' @rdname print
#' @method print COV
#' @aliases print.COV
#' @export
#'
print.COV <- function(x, ..., n=2L){
    stopifnot(length(n)==1L)
    if (is.array(x)) {
        n <- min(n, dim(x)[3])
        print(x[, , seq_len(n)])
        cat(" ... ")
        print(x[, , dim(x)[3] - seq_len(n)])
    } else {
        n <- min(n, length(x))
        print(utils::head(x, n))
        cat("\n...\n\n")
        print(utils::tail(x, n))
    }
}
#' 
#' @rdname print
#' @method print lrt
#' @aliases print.lrt
#' @export
#'
print.lrt <- function(x, ..., dist=c("n", "c")){
    dist <- match.arg(dist)
    x1 <- data.table::copy(x)
    data.table::setattr(x1, "class", "data.frame")
    rownames(x1) <- x1[, "W"]
    x1[, "W"] <- NULL
    if (ncol(x1)==3) {
        x1[, "pChisq"] <- format.pval(x1[, "pChisq"])
        stats::printCoefmat(x1,
                            has.Pvalue=TRUE,
                            cs.ind=1L, # *c*oefficients and *s*tandard errors
                            dig.tst=getOption("digits"))
    } else { 
        x1[, c("pNorm", "pChisq")] <- format.pval(x1[, c("pNorm", "pChisq")])
        ## no need to print pChiSq values routinely
        if (dist=="n") {
            stats::printCoefmat(x1[, c("Q", "Var", "Z", "pNorm")],
                                has.Pvalue=TRUE,
                                cs.ind=seq.int(2), 
                                dig.tst=getOption("digits"))
        } else {
            stats::printCoefmat(x1[, c("Q", "Var", "chiSq", "df", "pChisq")],
                                has.Pvalue=TRUE,
                                cs.ind=as.integer(c(1, 2)),
                                dig.tst=getOption("digits"))
        }
    }
}
#'
#' @rdname print
#' @method print sup
#' @aliases print.sup
#' @export
#'
print.sup <- function(x, ...){
    x1 <- data.table::copy(x)
    data.table::setattr(x1, "class", "data.frame")
    rownames(x1) <- x1[, "W"]
    x1[, "W"] <- NULL
    x1[, "pSupBr"] <- format.pval(x1[, "pSupBr"])
    stats::printCoefmat(x1,
                        has.Pvalue=TRUE,
                        cs.ind=seq.int(2), # *c*oefficients and *s*tandard errors
                        dig.tst=getOption("digits"))
}
#'
#' @rdname print
#' @method print tableAndPlot
#' @aliases print.tableAndPlot
#' @export
#'
print.tableAndPlot <- function(x, ...,
                               hideTabLeg=TRUE,
                               tabHeight=0.25){
    autoplot(x,
             hideTabLeg=hideTabLeg,
             tabHeight=tabHeight)
}
#'
#' @rdname print
#' @method print stratTableAndPlot
#' @aliases print.stratTableAndPlot
#' @export
#'
print.stratTableAndPlot <- function(x, ...,
                                    hideTabLeg=TRUE,
                                    tabHeight=0.25) {
    for (i in seq.int(length(x))) {
        if (interactive()) {
            ## max. no devices is 63
            if (i %% 63 == 0) grDevices::graphics.off()
            grDevices::dev.new()
        }
        autoplot(x[[i]],
                 hideTabLeg=hideTabLeg,
                 tabHeight=tabHeight)
    }
}
