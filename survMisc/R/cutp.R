#' @name cutp
#' @title \bold{cut p}oint for a continuous variable in a
#' model fit with \code{coxph} or \code{survfit}.
#'
#' @include ten.R
#' @include print.R
#' 
#' @description Determine the optimal cut point for a continuous variable
#' in a \code{coxph} or \code{survfit} model.
#'
#' @param x A \code{survfit} or \code{coxph} object
#' @param defCont \bold{def}inition of a \bold{cont}inuous variable.
#' \cr
#' If the variable has \eqn{>} \code{defCont} unique values, it
#' is treated as continuous and a cut point is determined. 
#' @param ... Additional arguments (not implemented).
#' 
#' @return A \code{list} of \code{data.table}s.
#'  \cr
#' There is one list element per continuous variable.
#'  \cr
#' Each has a column with possible values of the cut point
#'  (i.e. unique values of the variable), and the
#' additional columns:
#'  \item{U}{The score (log-rank) test for a model with the variable 'cut'
#'   into into those \eqn{\geq}{>=} the cutpoint and those below.}
#'  \item{Q}{The test statistic.}
#'  \item{p}{The \eqn{p}-value.}
#' The tables are ordered by \eqn{p}-value, lowest first.
#' 
#' @details
#' For a cut point \eqn{\mu}{mu}, of a predictor \eqn{K},
#' the variable is split
#' into two groups, those \eqn{\geq \mu}{>= mu} and
#' those \eqn{< \mu}{< mu}.
#'  \cr
#' The score (or log-rank) statistic, \eqn{sc},
#' is calculated for each unique element
#' \eqn{k} in \eqn{K} and uses
#' \itemize{
#'  \item \eqn{e_i^+}{e1[i]} the number of events
#'  \item \eqn{n_i^+}{n1[i]} the number at risk
#' }
#' in those above the cut point, respectively.
#'  \cr
#' The basic statistic is 
#' \deqn{sc_k = \sum_{i=1}^D ( e_i^+ - n_i^+ \frac{e_i}{n_i} )}{
#'       sc[k] = sum (e1[i] - n1[i] * e[i] / n[i])}
#'  \cr
#' The sum is taken across times with observed events, to \eqn{D},
#' the largest of these.
#'  \cr
#' It is normalized (standardized), in the case of censoring,
#' by finding \eqn{\sigma^2}{s^2} which is:
#' \deqn{\sigma^2 = \frac{1}{D - 1}
#'                  \sum_i^D (1 - \sum_{j=1}^i \frac{1}{D+ 1 - j})^2}{
#'       s^2 = (1 / (D - 1)) *
#'             sum[i:D](1 - sum[j:i](1 / (D - j + 1))^2 )}
#' The test statistic is then
#' \deqn{Q = \frac{\max |sc_k|}{\sigma \sqrt{D-1}}}{
#'       Q = max(abs(sc[k])) / s * sqrt((D - 1))}
#' Under the null hypothesis that the chosen cut point
#' does \emph{not} predict survival,
#' the distribution of \eqn{Q} has a limiting distibution which
#' is the supremum of the
#' absolute value of a Brownian bridge:
#' \deqn{p = Pr(\sup Q \geq q) = 2 \sum_{i=1}^{\infty}
#'                              (-1)^{i + 1} \exp (-2 i^2 q^2)}{
#'       p= P(Q >= q) = 2 sum[i:Inf](-1)^(i + 1) * e^(-2 * i^2 *q^2)}
#'
#' @references Contal C, O'Quigley J, 1999.
#' An application of changepoint methods in studying the
#' effect of age on survival in breast cancer.
#' \emph{Computational Statistics & Data Analysis} \bold{30}(3):253--70.
#' \href{http://dx.doi.org/10.1016/S0167-9473(98)00096-6}{
#'       ScienceDirect (paywall)}
#' 
#' @references Mandrekar JN, Mandrekar, SJ, Cha SS, 2003.
#' Cutpoint Determination Methods in Survival Analysis using SAS.
#' \emph{Proceedings of the 28th SAS Users Group International Conference (SUGI)}. Paper 261-28.
#' \href{http://www2.sas.com/proceedings/sugi28/261-28.pdf}{
#'       SAS (free)}
#'
#' @rdname cutp
#' @export
#'
cutp <- function(x, ...) UseMethod("cutp")
#'
#' @rdname cutp
#' @method cutp coxph
#' @aliases cutp.coxph
#' @export
#' @examples
#' ## Mandrekar et al. above
#' data("bmt", package="KMsurv")
#' b1 <- bmt[bmt$group==1, ] # ALL patients
#' c1 <- coxph(Surv(t2, d3) ~ z1, data=b1) # z1=age
#' c1 <- cutp(c1)$z1
#' data.table::setorder(c1, "z1")
#' ## [] below is used to print data.table to console
#' c1[]
#' 
#' \dontrun{
#' ## compare to output from survival::coxph
#' matrix(
#'     unlist(
#'         lapply(26:30,
#'                function(i) c(i, summary(coxph(Surv(t2, d3) ~ z1 >= i, data=b1))$sctest))),
#'     ncol=5,
#'     dimnames=list(c("age", "score_test", "df", "p")))
#' cutp(coxph(Surv(t2, d3) ~ z1, data=bmt[bmt$group==2, ]))$z1[]
#' cutp(coxph(Surv(t2, d3) ~ z1, data=bmt[bmt$group==3, ]))[[1]][]
#' ## K&M. Example 8.3, pg 273-274.
#' data("kidtran", package="KMsurv")
#' k1 <- kidtran
#' ## patients who are male and black
#' k2 <- k1[k1$gender==1 & k1$race==2, ]
#' c2 <- coxph(Surv(time, delta) ~ age, data=k2)
#' print(cutp(c2))
#' ## check significance of computed value
#' summary(coxph(Surv(time, delta) ~ age >= 58, data=k2))
#' k3 <- k1[k1$gender==2 & k1$race==2, ]
#' c3 <- coxph(Surv(time, delta) ~ age, data=k3)
#' print(cutp(c3))
#' ## doesn't apply to binary variables e.g. gender
#' print(cutp(coxph(Surv(time, delta) ~ age + gender, data=k1)))
#' }
#' 
cutp.coxph <- function(x, ...,
                       defCont=3){
    stopifnot(inherits(x, "coxph"))
    d1 <- attr(x$terms, "dataClasses")[-1]
    d1 <- d1[d1=="numeric"]
    m1 <- data.frame(stats::model.matrix(x))
    res1 <- lapply(seq.int(d1), function(i) {
        ## variable
        var1 <- get(names(d1[i]), m1)
        ## unique values
        u1 <- unique(var1)
        if (length(u1) < defCont) return(NaN)
        ## convert to logical
        l1 <- lapply(u1, function(j) var1 >= j)
        names(l1) <- u1
        
        ## get tne; this is the longest step
        l1 <- lapply(l1, function(j) asWide(ten(x$y ~ j)))
        ## U = score test
        l1 <- lapply(l1, function(dt1) dt1[, sum(e_1 - e * n_1 / n)])
        res2 <- data.table::data.table(
            u1,
            "U"=abs(unlist(l1)))
        data.table::setnames(res2, old="u1", new=names(d1[i]))
        data.table::setorder(res2, -U)
        s2 <- findS2(sum(x$nevent))
        res2[, "Q" := U / (sqrt(s2) * sqrt(sum(x$nevent) - 1))]
        res2[, "p" := unlist(lapply(Q, findP))]
        return(res2)
    })
    names(res1) <- names(d1)
    return(res1)
}
### helper functions
findS2 <- function(D) {
    (1 / (D - 1)) *
        sum(
            sapply(1:D,
                   ## i in 1:D
                   function(i)
                       ## j in 1:i
                   (1 - sum(sapply(1:i,
                                   function(j)
                                       1 / (D + 1 - j))))^2))
}
###
## lim = limit (accuracy)
##  should be to Inf but generally 1e3 is enough
findP <- function(q, lim=1e3) {
    if (q < 0.2) return(1)
    2 * sum(sapply(seq.int(lim),
                   function(j) {(-1)^(j+1) * exp(-2 * j^2 * q^2)}))
}
#'
#' @rdname cutp
#' @method cutp survfit
#' @aliases cutp.survfit
#' @export
#' 
cutp.survfit <- function(x, ..., defCont=3){
    f1 <- deparse(x$call)
    f1 <- sub("survfit", "coxph", f1)
    c1 <- eval(parse(text=f1))
    cutp(c1, defCont=defCont)
}
## R CMD check
U <- Q <- e_1 <- n_1 <- NULL
## if (plot){
##            m1 <- paste0(
##              "Test statistic for cut points \n For variable ", var,
##              "\nLarger values indicate cut point more likely here")
##            setkey(res1, var)
##            res1[, graphics::plot(var, U,
##                                  xlab="Cut point",
##                                  ylab="Test statistic",
##                                  main=m1,
##                                  ...)]
##            res1[, graphics::lines(var, U, ...)]
## }
