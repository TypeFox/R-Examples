#' @name dx
#' @title Diagnostics for binomial regression
#'
#' @rdname dx
#' @export
#'
dx <- function(x, ...){
    UseMethod("dx")
}
#' @rdname dx
#' @aliases dx.glm
#' @method dx glm
#' @export
#'
#' @description
#'  Returns diagnostic measures for a binary regression
#'  model by covariate pattern
#' @param x A regression model with class \code{glm} and
#'  \code{x$family$family == "binomial"}.
#' @param ... Additional arguments
#'  which can be passed to:
#'  \cr
#'  ?stats::model.matrix
#'  \cr
#'  e.g. \code{contrasts.arg} which can be used for \code{factor} coding.
#' @param byCov Return values by \emph{covariate pattern},
#'  rather than by individual observation.
#' @return A \code{data.table}, with rows sorted by \eqn{\Delta \hat{\beta}_i}{dBhat}.
#' \cr \cr
#' If \code{byCov==TRUE}, there is one row per covariate pattern
#' with at least one observation.
#' \cr \cr
#' The initial columns give the predictor variables \eqn{1 \ldots p}{1 ... p}.
#' \cr
#' Subsequent columns are labelled as follows:
#'
#' \item{\eqn{\mathrm{y} \quad y_i}{y}}{
#'  The \emph{actual} number of observations with
#'  \eqn{y=1} in the model data.
#' }
#'
#' \item{\eqn{\mathrm{P} \quad P_i}{P}}{
#'  Probability of this covariate pattern.
#'  \cr
#'  This is given by the  inverse of the link function,
#'  \code{x$family$linkinv}. See:
#'  \cr
#'  ?stats::family
#' }
#'
#' \item{\eqn{\mathrm{n} \quad n_i}{n}}{
#'  Number of observations with these covariates.
#'  \cr
#'  If \code{byCov=FALSE} then this will be \eqn{=1} for all observations.
#' }
#'
#' \item{\eqn{\mathrm{yhat} \quad \hat{y}}{yhat}}{
#'  The \emph{predicted} number of observations having
#'  a response of \eqn{y=1}, according to the model.
#'  \cr
#'  This is:
#'   \deqn{\hat{y_i} = n_i P_i}{
#'         yhat[i] = n[i] * P[i]}
#' }
#'
#' \item{\eqn{\mathrm{h} \quad h_i}{h}}{
#'  Leverage, the diagonal of the \bold{h}at matrix used to
#'  generate the model:
#'   \deqn{H = \sqrt{V} X (X^T V X)^{-1} X^T \sqrt{V}}{
#'         H = V^0.5 X (X'VX)^-1 X'V^0.5}
#'  Here \eqn{^{-1}}{^-1} is the inverse and
#'  \eqn{^T}{'} is the transpose of a matrix.
#'  \cr
#'  \eqn{X} is the matrix of predictors, given by \code{stats::model.matrix}.
#'  \cr
#'  \eqn{V} is an \eqn{N \times N}{N * N} sparse matrix. All elements are
#'  \eqn{=0} except for the diagonal, which is:
#'    \deqn{v_{ii} = n_iP_i (1 - P_i)}{
#'          v[i][i] = n[i]P[i] * (1 - P[i])}
#'  Leverage \eqn{H} is also the estimated covariance matrix of
#'  \eqn{\hat{\beta}}{Bhat}.
#'  \cr
#'  Leverage is  measure of the influence of this
#'  covariate pattern on the model and is approximately
#'   \deqn{h_i \approx x_i - \bar{x} \quad \mathrm{for} \quad 0.1 < P_i < 0.9}{
#'         h[i] = x[i] - mean(x), 0.1 < P[i] < 0.9}
#'  That is, leverage is approximately equal to the distance of
#'  the covariate pattern \eqn{i} from the mean \eqn{\bar{x}}{mean(x)}.
#'  \cr
#'  For values of \eqn{p} which are large (\eqn{>0.9}) or
#'  small (\eqn{<0.1}) this relationship no longer holds.
#' }
#'
#' \item{\eqn{\mathrm{Pr} \quad Pr_i}{Pr}}{
#'  The Pearson residual, a measure of influence. This is:
#'   \deqn{Pr_i = \frac{y_i - \mu_y}{\sigma_y}}{
#'      Pr[i] = (y[i] - ybar) / SD[y]}
#'  where \eqn{\mu_y}{ybar} and \eqn{\sigma_y}{SD[y]} refer
#'  to the mean and standard deviation of a binomial distribution.
#'  \cr
#'  \eqn{\sigma^2_y = Var_y}{SD^2 = Var}, is the variance.
#'   \deqn{E(y=1) = \mu_y = \hat{y} = nP \quad \mathrm{and} \quad \sigma_y=\sqrt{nP(1 - P)}}{
#'         E(y=1) = ybar = yhat = nP and SE[y] = (nP(1-P))^0.5}
#'  Thus:
#'   \deqn{Pr_i = \frac{y_i - n_i P_i}{\sqrt{n_i P_i (1 - P_i)}}}{
#'         Pr[i] = (y[i] - n[i]P[i]) / (n[i]P[i](1 - P[i])^0.5)}
#' }
#'
#' \item{\eqn{\mathrm{dr} \quad dr_i}{dr}}{
#'  The deviance residual, a measure of influence:
#'   \deqn{dr_i = \mathrm{sign}(y_i - \hat{y}_i) \sqrt{d_i}}{
#'         dr[i] = sign(y[i] - yhat[i]) * d[i]^0.5}
#'  \eqn{d_i}{d[i]} is the contribution of observation \eqn{i} 
#'  to the model deviance.
#'  \cr
#'  The \eqn{\mathrm{sign}}{sign} above is:
#'   \itemize{
#'    \item \eqn{y_i > \hat{y}_i \quad \rightarrow \mathrm{sign}(i)=1}{
#'               y[i] > yhat --> sign(i) = 1}
#'    \item \eqn{y_i = \hat{y}_i \quad \rightarrow \mathrm{sign}(i)=0}{
#'               y[i] = yhat --> sign(i) = 0}
#'    \item \eqn{y_i < \hat{y}_i \quad \rightarrow \mathrm{sign}(i)=-1}{
#'               y[i] < yhat --> sign(i) = -1}
#'   }
#'  In logistic regression this is:
#'   \deqn{y_i = 1 \quad \rightarrow \quad dr_i = \sqrt{2 \log (1 + \exp(f(x))) - f(x)}}{
#'         y[i] = 1 --> dr[i] = (2 * log (1 + e^f(x) - f(x)))^0.5}
#'   \deqn{y_i = 0 \quad \rightarrow \quad dr_i = -\sqrt{2 \log (1 + \exp(f(x)))}}{
#'         y[i] = 0 --> dr[i] = (2 * log (1 + e^f(x)))}
#'  where \eqn{f(x)} is the linear function of the predictors
#'  \eqn{1 \ldots p}{1 ... p}:
#'   \deqn{f(x) = \hat{\beta_0} + \hat{\beta_1} x_{1i} + \ldots + \hat{\beta_p} x_{ip}}{
#'         f(x) = B[0] + B[1][i] * x[1][i] + ... + B[p][i] * x[p][i]}
#'  this is also:
#'   \deqn{dr_i = sign(y_i - \hat{y_i}) \sqrt{
#'                 2 (y_i \log{\frac{y_i}{\hat{y_i}}} +
#'                 (n_i - y_i) \log{\frac{n_i - y_i}{n_i(1-p_i)}} )}}{
#'         dr[i] = sign(y[i] - yhat[i])
#'                 [2 * (y[i] * log(y[i] / n[i] * p[i])) +
#'                 (n[i] - y[i]) * log((n[i] - y[i]) / (n[i] * (1 - p[i])))]^0.5}
#' To avoid the problem of division by zero:
#'  \deqn{y_i = 0 \quad \rightarrow \quad dr_i = - \sqrt{2n_i| \log{1 - P_i} | }}{
#'        y[i] = 0 --> dr[i] = (2 * n[i] * | log(1 - P[i]) |)^0.5}
#' Similarly to avoid \eqn{\log{\infty}}{log(Inf)}:
#'  \deqn{y_i = n_i \quad \rightarrow \quad dr_i = \sqrt{2n_i | \log{P_i} | }}{
#'        y[i] = n[i] --> dr[i] = (2 * n[i] * | log(P[i]) |)^0.5}
#' The above equations are used when calculating \eqn{dr_i}{dr[i]} by covariate group.
#' }
#'
#' \item{\eqn{\mathrm{sPr} \quad sPr_i}{sPr}}{
#'  The \bold{s}tandardized Pearson residual.
#'  \cr
#'  The residual is standardized by the leverage \eqn{h_i}{h[i]}:
#'   \deqn{sPr_i = \frac{Pr_i}{\sqrt{(1 - h_i)}}}{
#'         sPr[i] = Pr[i] / (1 - h[i])^0.5}
#' }
#'
#' \item{\eqn{\mathrm{sdr} \quad sdr_i}{sdr}}{
#'  The \bold{s}tandardized deviance residual.
#'  \cr
#'  The residual is standardized by the leverage, as above:
#'   \deqn{sdr_i = \frac{dr_i}{\sqrt{(1 - h_i)}}}{
#'         sdr[i] = dr[i] / (1 - h[i])^0.5}
#' }
#'
#' \item{\eqn{\mathrm{dChisq \quad \Delta P\chi^2_i}}{dChisq}}{
#'  The change in the Pearson chi-square statistic
#'  with observation \eqn{i} removed.
#'  Given by:
#'   \deqn{\Delta P\chi^2_i = sPr_i^2 = \frac{Pr_i^2}{1 - h_i}}{
#'         dChi^2 = sPr[i]^2 = Pr[i]^2 / (1 - h[i])}
#'  where \eqn{sPr_i}{sPr[i]} is the standardized Pearson residual,
#'  \eqn{Pr_i}{Pr[i]} is the Pearson residual and
#'  \eqn{h_i}{h[i]} is the leverage.
#'  \cr
#'  \eqn{\Delta P\chi^2_i}{dChisq[i]} should be \eqn{<4}
#'  if the observation has little influence on the model.
#' }
#'
#' \item{\eqn{\Delta D_i  \quad \mathrm{dDev}}{dDev}}{
#'  The change in the deviance statistic
#'  \eqn{D = \sum_{i=1}^n dr_i}{D = SUM dr[i]}
#'  with observation \eqn{i} excluded.
#'  \cr
#'  It is scaled by the leverage \eqn{h_i}{h[i]} as above:
#'    \deqn{\Delta D_i = sdr_i^2 = \frac{dr_i^2}{1 - h_i}}{
#'           dDev[i] = sdr[i]^2 = dr[i]^2 / (1 - h[i])}
#' }
#'
#' \item{\eqn{\Delta \hat{\beta}_i \quad \mathrm{dBhat}}{dBhat}}{
#'  The change in \eqn{\hat{\beta}}{Bhat}
#'  with observation \eqn{i} excluded.
#'  \cr
#'  This is scaled by the leverage as above:
#'    \deqn{\Delta \hat{\beta} = \frac{sPr_i^2 h_i}{1 - h_i}}{
#'          dBhat = h[i] * sPr[i]^2 / (1 - h[i])}
#'  where \eqn{sPr_i}{sPR[i]} is the standardized Pearson residual.
#'  \cr
#'  \eqn{\Delta \hat{\beta}_i}{dBhat[i]} should be \eqn{<1}
#'  if the observation has little influence on the model coefficients.
#' }
#'
#' @note By default, values for the statistics are calculated by
#' \emph{covariate pattern}.
#' Different values may be obtained if
#' calculated for each individual
#' obervation (e.g. rows in a \code{data.frame}).
#' \cr \cr
#' Generally, the values calculated by
#' covariate pattern are preferred,
#' particularly where the number of observations in a group is \eqn{>5}.
#' \cr
#' In this case Pearsons chi-squared and the deviance statistic
#' should follow a chi-squared distribution with \eqn{i - p} degrees of freedom.
#' 
#' @seealso \code{\link{plot.glm}}
#' 
#' @examples
#' ## H&L 2nd ed. Table 5.8. Page 182.
#' ## Pattern nos. 31, 477, 468
#' data(uis)
#' uis <- within(uis, {
#'     NDRGFP1 <- 10 / (NDRGTX + 1)
#'     NDRGFP2 <- NDRGFP1 * log((NDRGTX + 1) / 10)
#' })
#' (d1 <- dx(g1 <- glm(DFREE ~ AGE + NDRGFP1 + NDRGFP2 + IVHX +
#'                     RACE + TREAT + SITE +
#'                     AGE:NDRGFP1 + RACE:SITE,
#'                     family=binomial, data=uis)))
#' d1[519:521, ]
dx.glm <- function(x, ..., byCov=TRUE){
    stopifnot(inherits(x, "glm"))
    stopifnot(x$family$family=="binomial")
### for R CMD check
    y <- P <- n <- yhat <- NULL
    Pr <- dr <- h <- sPr <- sdr <- NULL
    dBhat <- dChisq <- dDev <- NULL
    c1 <- x$coefficients
    ## get model matrix
    res1 <- data.table::data.table(stats::model.matrix(x))
    ## individual observations
    res1[, "y" := x$y]
    res1[, "P" := x$fitted.values]
    res1[, "n" := rep(1L, nrow(res1))]
    res1[, "yhat" := n * P]
    res1[, "Pr" := stats::residuals.glm(x, "pearson")]
    res1[, "dr" := stats::residuals.glm(x, "deviance")]
    res1[, "h" := stats::hatvalues(x)]
    ## unique co-variate patterns
    if(byCov){
        res1[, "y" := sum(y), by=P]
        res1[, "n" := length(n), by=P]
        ## remove duplicates
        res1 <- res1[!duplicated(P), ]
        res1[, "yhat" := n * P]
        res1[, "Pr" := (y - n * P) / (sqrt(n * P * (1 - P)))]
        .getDr <- function(y, n, P){
            if(y==0){
                return(- sqrt(2 * n * abs(log(1 - P))))
            } else if(y==n){
                return(sqrt(2 * n * abs(log(P))))
            } else return( sign(y - (n * P)) *
                sqrt(2 * ((y * log(y / (n * P))) + (n - y) * log((n - y) / (n * (1-P))))))
        }
        res1[, "dr" := .getDr(y, n, P), by=P]
    }
    res1[, "sPr" := Pr / sqrt(1 - h)]
    res1[, "sdr" := dr / sqrt(1 - h)]
    res1[, "dChisq" := sPr^2]
    res1[, "dDev" := dr^2 / (1 - h)]
    res1[, "dBhat" := (sPr^2 * h ) / (1 - h)]
###
    data.table::setkey(res1, dBhat)
    class(res1) <- c("dxBinom", class(res1))
    return(res1)
}
