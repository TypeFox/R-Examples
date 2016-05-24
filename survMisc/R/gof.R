#' @name gof
#' @title \bold{g}oodness \bold{o}f \bold{f}it test for a \code{coxph} object
#' 
#' @rdname gof
#' @export
#' 
gof <- function(x, ...)  UseMethod("gof")
#'
#' @rdname gof
#' @method gof coxph
#' @aliases gof.coxph
#' @export
#' 
#' @param x An object of class \code{coxph}
#' @param ... Additional arguments (not implemented)
#' @param G  Number of \bold{g}roups into which to divide risk score.
#' If \code{G=NULL} (the default), uses closest integer to
#'  \deqn{G = \max(2, \quad \min(10, \quad \frac{ne}{40}))}{
#'        G = max(2, min(10, ne/40))}
#' where \eqn{ne} is the number of events overall.
#'
#' @return A \code{list} with elements:
#' \item{groups}{A \code{data.table} with one row per group \eqn{G}.
#' The columns are \describe{
#'   \item{n}{Number of observations}
#'   \item{e}{Number of events}
#'   \item{exp}{Number of events expected. This is
#'              \deqn{exp = \sum e_i - M_i}
#'              where \eqn{e_i} are the events and
#'              \eqn{M_i} are the martingale residuals
#'               for each observation \eqn{i}}
#'   \item{z}{\eqn{Z} score, calculated as
#'     \deqn{ Z = \frac{e - exp}{\sqrt{exp}}}{
#'            Z = (e - exp) / exp^0.5}
#'  }
#'   \item{p}{\eqn{p}-value for \eqn{Z}, which is
#'      \deqn{ p = 2. \code{pnorm}(-|z|)}{
#'             p = 2 * pnorm(-|z|)}
#'    where \code{pnorm} is the normal distribution function
#'  with mean \eqn{\mu =0}{0} and standard deviation \eqn{\sigma =1}{1}
#' and \eqn{|z|} is the absolute value.}
#' }}
#' \item{lrTest}{Likelihood-ratio test.
#' Tests the improvement in log-likelihood with addition
#' of an indicator variable with \eqn{G-1} groups.
#' This is done with \code{survival:::anova.coxph}.
#' The test is distributed as chi-square with \eqn{G-1} degrees of freedom}
#' 
#' @details
#' In order to verify the overall goodness of fit,
#' the risk score \eqn{r_i}{r[i]} for each observation \eqn{i} is given by
#' \deqn{r_i = \hat{\beta} X_i}{r[i] = B.X[i]}
#' where \eqn{\hat{\beta}}{B} is the vector of fitted coefficients
#' and \eqn{X_i}{X[i]} is the vector of predictor variables for
#' observation \eqn{i}.
#'  \cr
#' This risk score is then sorted and 'lumped' into
#' a grouping variable with \eqn{G} groups,
#' (containing approximately equal numbers of observations).
#'  \cr
#' The number of observed (\eqn{e}) and expected (\eqn{exp}) events in
#' each group are used to generate a \eqn{Z} statistic for each group,
#' which is assumed to follow a normal distribution with
#' \eqn{Z \sim N(0,1)}.
#'  \cr
#' The indicator variable \code{indicG} is added to the
#' original model and the two models are compared to determine the
#' improvement in fit via the likelihood ratio test.
#'
#' @note The choice of \eqn{G} is somewhat arbitrary but rarely should
#' be \eqn{> 10}.
#'  \cr
#' As illustrated in the example, a larger value for
#' \eqn{G} makes the \eqn{Z} test for each group more likely to be significant.
#' This does \emph{not} affect the significance of adding the
#' indicator variable \code{indicG} to the original model.
#'  \cr \cr
#' The \eqn{Z} score is chosen for simplicity, as for large sample sizes
#' the Poisson distribution approaches the normal. Strictly speaking,
#' the Poisson would be more appropriate for \eqn{e} and \eqn{exp}{exp} as
#' per Counting Theory.
#'  \cr
#' The \eqn{Z} score may be somewhat conservative as the expected events
#' are calculated using the martingale residuals from the overall model,
#' rather than by group. This is likely to bring the expected events
#' closer to the observed events.
#'  \cr \cr
#' This test is similar to the Hosmer-Lemeshow test for logistic regression.
#' 
#' @source
#' Method and example are from: \cr
#' May S, Hosmer DW 1998.
#' A simplified method of calculating an overall goodness-of-fit test
#' for the Cox proportional hazards model.
#' \emph{Lifetime Data Analysis} \bold{4}(2):109--20.
#' \href{http://dx.doi.org/10.1023/A:1009612305785}{Springer (paywall)}
#' 
#' @references
#' Default value for \eqn{G} as per: \cr
#' May S, Hosmer DW 2004.
#' A cautionary note on the use of the Gronnesby and Borgan
#' goodness-of-fit test for the Cox proportional hazards model.
#' \emph{Lifetime Data Analysis} \bold{10}(3):283--91.
#' \href{http://dx.doi.org/10.1023/B:LIDA.0000036393.29224.1d}{Springer (paywall)}
#' @references
#' Changes to the \code{pbc} dataset in the example are as detailed in: \cr
#' Fleming T, Harrington D 2005.
#' \emph{Counting Processes and Survival Analysis}.
#' New Jersey: Wiley and Sons. Chapter 4, section 4.6, pp 188.
#' \href{http://dx.doi.org/10.1002/9781118150672}{Wiley (paywall)}
#'
#' @examples
#' data("pbc", package="survival")
#' pbc <- pbc[!is.na(pbc$trt), ]
#' ## make corrections as per Fleming
#' pbc[pbc$id==253, "age"] <-  54.4
#' pbc[pbc$id==107, "protime"] <-  10.7
#' ### misspecified; should be log(bili) and log(protime) instead
#' c1 <- coxph(Surv(time, status==2) ~
#'             age + log(albumin) + bili + edema + protime,
#'             data=pbc)
#' gof(c1, G=10)
#' gof(c1)
#' 
gof.coxph <- function(x,
                      ...,
                      G=NULL){
    stopifnot(inherits(x, "coxph"))
    if(is.null(G)) G <- round(max(2, min(10, x$nevent/40)), 0)
    ## from survival:::predict.coxph
    r1 <- stats::predict(x, type="risk")
    r2 <- r1[order(r1)]
    ## 'chunk' size
    csize1 <- length(r1) / G
    s1 <- split(r2, ceiling(seq_along(r2) / csize1))
    ## 'quantiles' (approx. equal no.s in each)
    q1 <- rep(1:G, unlist(lapply(s1, length)))
    ## events
    e1 <- stats::model.frame(x)[, 1][, "status"]
    dt2 <- data.table::data.table(n=tapply(e1[order(r1)], q1, length))
    dt2[, e := tapply(e1[order(r1)], q1, sum)]
    ## from survival:::predict.coxph
    m1 <- stats::residuals(x, type="martingale")
    ## cumulative hazard = events - margtingale
    cumHaz <- e1[order(r1)] -  m1[order(r1)]
    ## no. expected
    dt2[, exp:= tapply(cumHaz, q1, sum)]
    ## z = eerved - expected/sqrt(expected)
    dt2[, z := (e - exp) / sqrt(exp)]
    ## z-score
    dt2[, p := 2 * pnorm(-abs(z))]
    ## reformulate
    f1 <- paste0(as.character(x$formula)[3], " + indicG")
    f2 <- stats::as.formula(paste(x$formula[2], x$formula[1], f1))
    environment(f2) <- environment(x$formula)
    d1 <- as.character(x$call$data)
    d1 <- data.table::data.table(get(d1))
    d1 <- d1[order(r1), ]
    d1[, indicG := as.factor(q1)]
    c1 <- survival::coxph(formula=f2, data=d1)
    a1 <- stats::anova(x, c1)
    ## result
    res1 <- vector(mode="list", length=2)
    res1[[1]] <- dt2
    res1[[2]] <- a1
    names(res1) <- c("groups", "lrTest")
    return(res1)
}
## for R CMD check
e <- z <- p <- indicG <- NULL
