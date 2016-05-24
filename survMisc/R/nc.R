#' @name nc
#' @title Add \bold{n}umber \bold{c}ensored.
#'
#' @include ten.R
#' @include print.R
#' 
#' @param x An object of class \code{ten} or \code{stratTen}.
#' @inheritParams sf.ten
#' 
#' @return
#' The original object, with new column(s) added indicating the 
#' number censored at each time point, depending on \code{attr(x, "shape")}: 
#' \item{"long"}{the new column, \code{c}, gives
#'  the number censored at each timepoint, by covariate group.}
#' \item{"wide"}{new columns, beginning with \code{c_}, give
#'  the number censored at each timepoint, by covariate group. 
#'  There is an additional \code{nc} column giving 
#'  the \emph{total} number censored at each timepoint.}
#' A \code{stratTen} object has each \code{ten} element in the 
#' \code{list} modified as above.
#' 
#' @rdname nc
#' @export
#' 
nc <- function(x, ...) UseMethod("nc")
#'
#' @rdname nc
#' @method nc ten
#' @aliases nc.ten
#' @export
#' @examples
#' data("kidney", package="KMsurv")
#' t1 <- ten(survfit(Surv(time, delta) ~ type, data=kidney))
#' nc(t1)
#' nc(asWide(t1))
#' 
nc.ten <- function(x, ...){
    if (attr(x, "shape")=="long") {
        x[, "nc" := (c(-diff(ncg), tail(ncg, 1))-e), by=cg]
    } else {
        n_ <- grep("n_", names(x))
        e_ <- grep("e_", names(x))
        ## no. at risk - no. events
        nMe1 <- x[, .SD, .SDcols=n_] - x[, .SD, .SDcols=e_]
        ## add no. censored
        c1 <- nMe1[seq.int(nrow(x) - 1), ] -
            x[seq.int(2, nrow(x)), .SD, .SDcols=n_]
        c1 <- data.table::rbindlist(list(
            c1,
            x[nrow(x), .SD, .SDcols=n_]))
        ## names for censored columns
        c_ <- grep("e_", names(x), value=TRUE)
        substr(c_, 1, 1) <- "c"
        x[, (c_) := c1]
        ## total no. censored per time period
        x[, "c" := rowSums(.SD),
          .SDcols = grep("c_", names(x))]
        ## reorder
        data.table::setcolorder(x,
                                c("t", "n", "e", "c", colnames(x)[4L:(ncol(x) - 1L)]))
    }
    return(x)
}
#' @rdname nc
#' @method nc stratTen
#' @aliases nc.stratTen
#' @export
#' @examples
#' ## stratified model
#' data("pbc", package="survival")
#' t1 <- ten(coxph(Surv(time, status==2) ~ log(bili) + age + strata(edema), data=pbc))
#' nc(t1)
#' 
nc.stratTen <- function(x, ...){
    lapply(x, nc)
    return(x)
}
