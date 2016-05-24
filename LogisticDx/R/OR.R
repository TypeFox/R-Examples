#' @name OR
#' @title \bold{O}dds \bold{r}atio
#' for binary regression models fit with \code{glm}
#'
#' @rdname OR
#' @export
#'
OR <- function(x, ...){
    UseMethod("OR")
}
#' @rdname OR
#' @aliases OR.default
#' @method OR default
#' @export
OR.default <- function(x, ...){
    stopifnot(x >= 0 && x <= 1)
    return(x / (1-x))
}
#'
#' @param x A numeric object
#' containing probabilities \eqn{P}.
#' \cr
#' I.e. the range of \eqn{P} must be \eqn{0} to \eqn{1}.
#' \cr
#' The odds ratio \eqn{OR} is given by:
#' \deqn{OR_i = \frac{P_i}{1-P_i} = \frac{\frac{P_1}{1-P_1}}{\frac{P_0}{1-P_0}} =
#'              \frac{\mathrm{odds_1}}{\mathrm{odds_0}}}{
#'       OR = P[i] / (1-P[i]) = odds[1] / odds[0]}
#'
#' There is a method for regression models with \code{class(x)==glm} and
#' \code{x$family$family == "binomial"}.
#' @param ... Not used.
#' @param newdata A \code{vector} of new variables to use.
#' \cr
#' There should be one value, in sequence, for each coefficient
#' in the model.
#' \cr
#' By default, values are calculated for a change in the value of the
#' coefficient for the predictor from \eqn{0} to \eqn{1}.
#' \cr
#' For continuous predictors changes of \eqn{> 1}
#' unit may have more practical significance.
#' @param ci If \code{ci=TRUE} (the default), include
#' a confidence interval for \eqn{P_i}{P[i]} and \eqn{OR_i}{OR[i]} in the returned values.
#' @param alpha Used to cacluate the confidence interval, which is:
#'  \deqn{\mathrm{CI} = x \pm Z_{1 - \alpha} \sigma}{
#'        CI = x +- Z[1-alpha].SD}
#' where the normal distribution \eqn{Z \sim N (0,1)}{Z = N(0,1)} and
#' \eqn{\sigma}{SD} is the standard deviation.
#' @param what See \bold{Value} below.
#'
#' @return A \code{data.table}.
#' Columns give the model, the value of the link function and
#' the associated probability \eqn{P_i}{P[i]} and odds ratio \eqn{OR_i}{OR[i]}.
#' \cr \cr
#' If \code{ci=TRUE}, will also give upper and lower bounds of the
#' confidence intervals for these values.
#' \cr \cr
#' Rows are determined by \code{what}:
#' \item{\code{what="model"}}{
#'  The value of the link function is given
#'  for the full model.
#'  \cr
#'  If an intercept term is included, the value if given
#'  with \emph{and} without the intercept.
#' }
#'
#' \item{\code{what="all"}}{
#'  The value of the link function is given for each
#'  \emph{combination} of coefficients in the model.
#' }
#'
#' \item{\code{what="data"}}{
#'  The value of the link function is given
#'  for each set of predictors in the data with which
#'  the model was fit.
#'  \cr
#'  This option will ignore the argument \code{newdata}.
#' }
#' @note
#' In the model formulas, the intercept term is specified as
#' \code{0} (absent) or \code{1} (present).
#' \cr
#' The variance of the values of the link function is:
#' \deqn{\sigma ^2 = \sum x_i^2 \sigma^2(\hat{\beta_i}) +
#'                   \sum 2x_ix_j cov(\hat{\beta_i}, \hat{\beta_j})}{
#'       Var = SUM(x[i]^2.var(B[i])) + SUM(2.x[i][x[j]].cov(B[i], B[j]))}
#' where \eqn{\sigma^2}{Var} is the variance and
#' \eqn{cov} is the covariance.
#' \cr \cr
#' @seealso
#' ?stats::predict.glm
#'
#' @rdname OR
#' @aliases OR.glm
#' @method OR glm
#' @export
#'
#' @examples
#' if(require("graphics")){
#'     plot(x <- seq(from=0.1, to=0.9, by=0.05), y=OR(x))}
#' ## H&L 2nd ed. Table 1.3. Page 10.
#' data(ageChd)
#' summary(g1 <- glm(chd ~ age, data=ageChd, family=binomial))
#' OR(g1)
#' attributes(OR(g1))
#' ## Table 1.4. Page 20.
#' stats::vcov(g1)
#' ## Table 2.3. Page 38.
#' data(lbw)
#' summary(g1 <- glm(LOW ~ LWT + RACE, data=lbw, family=binomial))
#' ## Table 2.4. Page 42.
#' vcov(g1)
#' ageChd$gr54 <- ageChd$age > 54
#' OR(glm(chd ~ gr54, data=ageChd, family=binomial))
#'
OR.glm <- function(x, ...,
                   newdata=rep(1L, length(stats::coef(x))),
                   ci=TRUE,
                   alpha=0.95,
                   what=c("model", "all", "data")){
    stopifnot(inherits(x, "glm"))
    stopifnot(x$family$family %in% c("binomial", "quasi"))
    stopifnot(alpha <= 1 && alpha >=0)
    stopifnot(length(newdata) == length(stats::coef(x)))
    ## for R CMD check
    link <- low <- up <- P <- lowP <- upP <- NULL
    ## coefficients and values (multipliers)
    c1 <- x$coefficients
    v1 <- c1 * newdata
    ## rename intercept
    if ((i1 <- "(Intercept)" %in% names(c1))) names(c1)[1] <- "1"
    what <- match.arg(what)
###
    if(what=="model"){
        link1 <- sum(v1)
        if (i1) link1 <- c(sum(v1[-1], na.rm=FALSE), link1)
        n1 <- paste(names(c1), collapse="+")
        if (i1) {
            names(v1)[1] <- "0"
            n1 <- c(paste(names(v1), collapse="+"), n1)
        }
        res1 <- data.table::data.table(n1, link1)
    }
    if(what=="all"){
        ## combination Function nCk
        ## binomial expansion; gives no. of values e.g. logits
        ## can also use ch1 <- 1:sum(choose(n=length(v1), k=seq_along(v1)))
        ch1 <- seq_along(v1)
        res1 <- data.table::data.table(
            unlist(
                sapply(ch1,
                       function(m) apply(utils::combn(names(v1), m), 2,
                                         FUN=paste, collapse="+"),
                       simplify=FALSE)),
            unlist(
                sapply(ch1,
                       function(m) colSums(utils::combn(v1, m)),
                       simplify=FALSE))
            )
    }
    if(what=="data"){
        mm1 <- stats::model.matrix(x)
        colnames(mm1) <- names(c1)
        f1 <- function(r) paste(colnames(mm1), r, sep="=", collapse="+")
        n1 <- unname(apply(mm1, 1, f1))
        link1 <- unname(x$linear.predictors)
        res1 <- data.table::data.table(n1, link1)
    }
    data.table::setnames(res1, c("model", "link"))
    if(ci){
        ci1 <- - (alpha - 1) / 2
        z1 <- abs(stats::qnorm(ci1))
        vcov1 <- stats::vcov(x)
        ## variance
        var1 <- sum(diag(vcov1) * newdata^2
                    ) + sum(
                        2 * vcov1[lower.tri(vcov1)] *
                        apply(utils::combn(newdata, 2), MARGIN=2, FUN=prod))
        se1 <- sqrt(var1)
        err1 <- se1 * z1
        res1[, "low" := link - err1]
        res1[, "up" := link + err1]
    }
    ## function = inverse of link;
    ## gives estimate of probability
    res1 <- res1[, "P" := x$family$linkinv(link)]
    if(ci){
        res1[, "lowP" := x$family$linkinv(low)]
        res1[, "upP" := x$family$linkinv(up)]
    }
    res1[, "OR" := P / (1 - P)]
    if(ci){
        res1[, "lowOR" := lowP / (1 - lowP)]
        res1[, "upOR" := upP / (1 - upP)]
        res1[, c("low", "up") := NULL]
    }
    ## change to name of link
    data.table::setnames(res1, c(names(res1)[1],
                                 x$family$link,
                                 names(res1)[3:length(res1)]))
    class(res1) <- c("OR", class(res1))
    data.table::setattr(res1, "what", what)
    data.table::setattr(res1, "ci", ci)
    if(ci) data.table::setattr(res1, "alpha", alpha)
    return(res1)
}
