#' @name predict
#' @title predicted events
#' 
#' @include ten.R
#' @include print.R
#' @include asWide.R
#'
#' @param object An object of class \code{ten}.
#' @param eMP Add column(s) indicating
#'  \bold{e}vents \bold{m}inus \bold{p}redicted.
#' @inheritParams sf.ten
#'
#' @return An \code{attribute}, \code{pred} is added
#' to \code{object}:
#'  \item{t}{Times with at least one observation}
#'  \item{P_}{\bold{p}redicted number of events}
#' And if \code{eMP==TRUE} (the default):
#'  \item{eMP_}{\bold{e}vents \bold{m}inus \bold{p}redicted}
#' The names of the \code{object}'s covariate groups are
#' used to make the suffixes of the column names (i.e. after the
#' \code{_} character).
#' 
#' @details
#' With \eqn{K} covariate groups, We use \eqn{ncg_{ik}}{ncg[i, k]},
#' the number at risk for group \eqn{k},
#' to calculate the number of expected events:
#'  \deqn{P_{ik} = \frac{e_i(ncg_{ik})}{n_i} \quad k=1, 2 \ldots K}{
#'        P[i, k] = e[i] * ncg[i, k] / n[i]}
#' 
#' @note There is a predicted value for each unique time, for each covariate group.
#'
#' @seealso
#' ?survival::predict.coxph
#' methods("predict")
#'
#' @rdname predict
#' @method predict ten
#' @aliases predict.ten
#' @export 
#' @examples
#' ## K&M. Example 7.2, Table 7.2, pp 209-210.
#' data("kidney", package="KMsurv")
#' k1 <- ten(Surv(time=time, event=delta) ~ type, data=kidney)
#' predict(k1)
#' predict(asWide(k1))
#' stopifnot(predict(asWide(k1))[, sum(eMP_1 + eMP_2)] <=
#'           .Machine$double.neg.eps)
#' ## Three covariate groups
#' ## K&M. Example 7.4, pp 212-214.
#' data("bmt", package="KMsurv")
#' b1 <- ten(Surv(time=t2, event=d3) ~ group, data=bmt)
#' predict(b1)
#' ## one group only
#' predict(ten(Surv(time=t2, event=d3) ~ 1, data=bmt))
#' 
predict.ten <- function(object, ...,
                        eMP=TRUE,
                        reCalc=FALSE){
    if (!reCalc & !is.null(attr(object, "pred"))) return (attr(object, "pred"))
    stopifnot(attr(object, "ncg")>=1)
    if (attr(object, "shape")=="long" & (attr(object, "ncg")==1)) {
        res1 <- object[, "P" := e * ncg / n][, list(t, P)]
        if(eMP) res1[, "eMP" := object[, e] - P]
    } else {
        x1 <- if(attr(object, "shape")=="wide") object else asWide(object)
        ## names of columns
        n1 <- x1[, as.matrix(.SD), .SDcols=grep("n_", names(x1))]
        res1 <- data.table::data.table(n1 * x1[, e / n])
        data.table::setnames(res1, paste("P", seq.int(ncol(res1)), sep="_"))
        if (eMP) {
            e1 <- x1[, as.matrix(.SD), .SDcols=grep("e_", names(x1))]
            na1 <- paste("eMP", seq.int(ncol(res1)), sep="_")
            res1[, (na1) := data.frame(as.matrix(res1) - e1)]
        }
    }
    data.table::setattr(res1, "class", c("pred", class(res1)))
    data.table::setattr(object, "pred", res1)
    return(attr(object, "pred"))
}
