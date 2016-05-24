#' @name sig
#' @title Significance tests for a
#' binary regression models fit with \code{glm}
#'
#' @rdname sig
#' @export
sig <- function(x, ...){
    UseMethod("sig")
}
#' @rdname sig
#' @aliases sig.glm
#' @method sig glm
#' @export
#'
#' @param x A regression model with class \code{glm} and
#' \code{x$family$family == "binomial"}.
#' @param ... Not used.
#' @param test What to test.
#' \itemize{
#'  \item If \code{test="var"} (the default),
#'        will test significance for
#'        each \emph{variable} in the model.
#'        \cr
#'        This includes the intercept, if present.
#'        \cr
#'        This means \code{factor}s are tested for
#'        \emph{all} \code{level}s simultaneously.
#'  \item If \code{test="coef"}, will test
#'        significance for each \emph{coefficient} in the model.
#'        \cr
#'        This means the 'dummy variables' created from
#'        \code{factor}s will be tested individually.
#' }
#'
#' @return
#' A \code{list} of \code{data.table}s as follows:
#'
#' \item{Wald}{The Wald test for each coefficient which is:
#'   \deqn{W = \frac{\hat{\beta}}{\hat{SE_{\beta}}}}{
#'         W = B / SE[B]}
#'   This should be normally distributed.}
#' \item{LR}{The \bold{l}ikelihood \bold{r}atio test
#'   for each coefficient:
#'   \deqn{LR = -2 \log{\frac{\mathrm{likelihood \phantom{+} without \phantom{+} variable}}{
#'                            \mathrm{likelihood \phantom{+} with \phantom{+} variable}}}}{
#'         LR = -2 * log(likelihood without / likelihood with variable)}
#'   which is:
#'   \deqn{LR = -2 \sum_{i=1}^n(y_i \log{\frac{P_i}{y_i}} +
#'            (1 - y_i) \log{\frac{1 - P_i}{1-y_i}})}{
#'         LR = -2 * SUM(y * log(P / y) +
#'            (1 - y) * log((1 - P) / (1 - y)))}
#' When comparing a fitted model to a saturated model
#' (i.e. \eqn{P_i = y_i}{P[i]=y[i]} and likelihood \eqn{=1}),
#' the \eqn{LR} is referred to as the model \emph{deviance},
#' \eqn{D}.
#' }
#'
#' \item{score}{The score test, also known as the
#'   Rao, Cochran-Armitage trend and the Lagrange multiplier test.
#'   \cr
#'   This removes a variable from the model, then assesses
#'   the change. For logistic regression this is based on:
#'    \deqn{\bar{y} = \frac{\sum_{i=1}^n y_i}{n}}{
#'          ybar = (SUM y[i]) / n}
#'   and
#'    \deqn{\bar{x} = \frac{\sum_{i=1}^n x_i n_i}{n}}{
#'          xbar = (SUM x[i] * n[i]) / n}
#'   The statistic is:
#'    \deqn{\mathrm{ST} = \frac{\sum_{i=1}^n y_i(x_i - \bar{x})}{
#'                           \sqrt{\bar{y}(1-\bar{y})\sum_{i=1}^n (x_i-\bar{x}^2)}}}{
#'          ST = SUM x[i](y[i] - ybar) / (ybar(1 - ybar) SUM (x[i] - xbar)^2)^0.5}
#'  If the value of the coefficient is correct, the test should follow
#'  a standard normal distribution.
#' }
#'
#'
#' @note
#' The result has the \code{class} \code{"sig.glm"}.
#' The \code{print} method for this \code{class} shows only
#' the model coefficients and \eqn{p} values.
#'
#' @seealso
#' ?aod::wald.test
#' \cr
#' ?statmod::glm.scoretest
#' \cr
#' For corrected score tests:
#' \cr
#' ?mdscore::mdscore
#'
#' @examples
#' data(ageChd)
#' ## H&L 2nd ed. Table 1.3. Page 10.
#' summary(g1 <- glm(chd ~ age, data=ageChd, family=binomial))
#' sig(g1)
#' data(lbw)
#' ## Table 2.2. Page 36.
#' summary(g2 <- glm(LOW ~ AGE + LWT + RACE + FTV,
#'                   data=lbw, family=binomial))
#' sig(g2)
#' ## Table 2.3. Pages 38-39.
#' summary(g3 <- glm(LOW ~ LWT + RACE,
#'                   data=lbw, family=binomial))
#' sig(g3, test="coef")
#' ## RACE is more significant when dropped as a factor
#' ## 
#' sig(g3, test="var")
sig.glm <- function(x, ..., test=c("var", "coef")){
    stopifnot(inherits(x, "glm"))
    stopifnot(x$family$family=="binomial")
### for R CMD check
    chiSq <- df <- Z <- pVal <- NULL
    test <- match.arg(test)
    res1 <- vector(mode="list")
###
    if (test=="coef"){
        mm1 <- stats::model.matrix(x)
    } else {
        if ("(Intercept)" %in% names(stats::coef(x))){
            mm1 <- cbind("(Intercept)"=1,
                         stats::model.frame(x)[, -1, drop=FALSE])
        } else {
            mm1 <- stats::model.frame(x)[, -1]
        }
    }
    if(!is.data.frame(mm1)){
        n1 <- colnames(mm1)
        mm1 <- data.frame(mm1)
        colnames(mm1) <- n1
    }
### Wald
    res1$Wald <- data.table("coef"=colnames(mm1),
                            "stat"=rep("Z", ncol(mm1)),
                            "val"=rep(0.1, ncol(mm1)),
                            "df"=rep(0.1, ncol(mm1)),
                            "pVal"=rep(0.1, ncol(mm1)))
    s1 <- summary(x)$coefficients[, c("z value", "Pr(>|z|)")]
    for(i in seq.int(ncol(mm1))){
        if(!is.na(m1 <- match(colnames(mm1)[i], rownames(s1)))){
            set(res1$Wald, i=i, j=3L, value=s1[m1, "z value"])
            set(res1$Wald, i=i, j=4L, value=NA)
            set(res1$Wald, i=i, j=5L, value=s1[m1, "Pr(>|z|)"])
        } else {
            g1 <- grep(paste0("^", colnames(mm1)[i]),
                       rownames(s1))
            a1 <- aod::wald.test(b=stats::coef(x),
                                 Sigma=stats::vcov(x),
                                 Terms=g1)$result$chi2
            set(res1$Wald, i=i, j=2L, value="chiSq")
            set(res1$Wald, i=i, j=3L, value=a1["chi2"])
            set(res1$Wald, i=i, j=4L, value=as.integer(a1["df"]))
            set(res1$Wald, i=i, j=5L, value=a1["P"])
        }
    }
### likelihood ratio
    res1$LR <- data.table(coef=colnames(mm1),
                          "chiSq"=(a1 <- rep(0.1, ncol(mm1))),
                          "df"=a1,
                          "pVal"=a1)
    ## log likehihood for full model
    ll1 <- exp(x$deviance / -2)
    ## ## degrees of freedom for null model
    ## dfn1 <- x$df.null
    for (i in seq.int(ncol(mm1))){
        m1 <- stats::glm(x$y ~ 0 + ., data=mm1[, -i, drop=FALSE],
                         family=stats::binomial)
        ll2 <- exp(m1$deviance / -2)
        set(res1$LR, i=i, j=1L, value=colnames(mm1)[i])
        set(res1$LR, i=i, j=2L, value=2 * log(ll1 / ll2))
        set(res1$LR, i=i, j=3L,
            value=ifelse(is.factor(mm1[, i]),
            length(levels(mm1[, i])) - 1,
            1))
    }
    res1$LR[, "pVal" := 1 - stats::pchisq(chiSq, df)]
### score
    res1$score <- data.table(coef=res1$LR$coef,
                             "Z"=a1,
                             "pVal"=a1)
    for (i in 1:ncol(mm1)){
        m2 <- stats::glm(x$y ~ 0 + ., data=mm1[, -i, drop=FALSE],
                         family=stats::binomial)
        if (is.factor(mm1[, i])){
            z1 <- statmod::glm.scoretest(
                m2,
                stats::model.matrix( ~ 0 + mm1[, i]))
            z1 <- sum(abs(z1))
        } else {
            z1 <- statmod::glm.scoretest(m2, mm1[, i])
        }
        set(res1$score, i=i, j=2L, value=z1)
    }
    res1$score[, "pVal" := stats::pnorm(abs(Z),
                      lower.tail=FALSE) * 2]
    class(res1) <- c("sig.glm", class(res1))
    attr(res1, "test") <- test
    attr(res1, "coef") <- stats::coef(x)
    return(res1)
}
