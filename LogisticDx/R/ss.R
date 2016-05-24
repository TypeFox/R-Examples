#' @name ss
#' @title \bold{S}ample \bold{s}ize for a given coefficient and events per
#' covariate for model
#'
#' @rdname ss
#' @export
ss <- function(x, ...){
    UseMethod("ss")
}
#' @rdname ss
#' @aliases ss.glm
#' @method ss glm
#' @export
#'
#' @param x A regression model with class \code{glm} and
#' \code{x$family$family == "binomial"}.
#' @param ... Not used.
#' @param alpha significance level \eqn{\alpha}{alpha}
#' for the null-hypothesis significance test.
#' @param beta power \eqn{\beta}{Beta} for the null-hypothesis significance test.
#' @param coeff Name of coefficient (variable) in the model to be tested.
#' @param std Standardize the coefficient?
#'  \cr
#'  If \code{std=TRUE} (the default), a continuous coefficent will be
#'  standardized, using the mean \eqn{\bar{x}}{xbar} and standard deviation
#'  \eqn{\sigma_x}{SD[x]}:
#'   \deqn{z_x = \frac{x_i - \bar{x}}{\sigma_x}}{
#'         z[x] = (x[i] - xbar) / SD[x]}
#' @param alternative The default, \code{alternative="one.sided"},
#'  checks the null hypothesis with \code{z = 1 - alpha}.
#'  \cr
#'  If \code{alternative="two.sided"}, \code{z = 1 - alpha/2} is used instead.
#' @param OR Odds ratio. The size of the change in the probability.
#' @param Px0 The probability that \eqn{x=0}.
#'  \cr
#'  If not supplied, this is estimated from the data.
#' 
#' @return A list of:
#'  \item{ss}{Sample size required to show coefficient for
#'   predictor is as given in the model rather than the alternative (by default \eqn{=0}).
#' }
#'  \item{epc}{Events per covariate; should be \eqn{>10} to make meaningful
#'   statements about the coefficients obtained.
#' }
#' 
#' @details
#' Gives the sample size necessary to demonstrate that a coefficient
#' in the model for the 
#' given predictor is equal to its given value
#' rather than equal to zero (or, if \code{OR} is supplied,
#' the sample size needed to check for such a change in probability).
#' \cr \cr
#' Also, the number of events per predictor.
#' \cr
#' This is the \emph{smaller} value of the outcome \eqn{y=0} and outcome \eqn{y=1}. 
#' \cr \cr
#' For a \bold{continuous} coefficient, the calculation uses
#' \eqn{\hat{\beta}}{Bhat}, the estimated coefficient from the model, 
#' \eqn{\delta}{delta}:
#'  \deqn{\delta = \frac{1 + (1 + \hat{\beta}^2) \exp{1.25\hat{\beta}^2}}{
#'                       1 + \exp{-0.25 \hat{\beta}^2}}}{
#'        delta = (1 + (1 + Bhat^2)exp(1.25 * Bhat^2)) / (1 + exp(1 + exp(-0.25 * Bhat^2)))}
#' and \eqn{P_0}{P[0]}, the probability calculated from the intercept term
#' \eqn{\beta_0}{B[0]} from the logistic model
#' \cr
#' \code{glm(x$y ~ coeff, family=binomial)}
#' \cr as
#'  \eqn{P_0 = \frac{\exp{\beta_0}}{1 + \exp{\beta_0}}}{
#'       P[0] = exp(B[0]) / (1 + exp(B[0]))}
#' For a model with one predictor, the calculation is:
#'  \deqn{n = (1 + 2P_0 \delta)
#'            \frac{z_{1-\alpha} + z_{\code{beta}} \exp{0.25 \hat{\beta}^2}^2}{
#'                  P_0 \hat{\beta}^2}}{
#'        n = (1 + 1 * P[0] * delta) *
#'        (z[1-alpha] + z[beta] exp((0.25 * Bhat)^2)^2) /
#'        P[0] * Bhat^2}
#' For a multivariable model, the value is adjusted by \eqn{R^2}, the correlation
#' of \code{coeff} with the other predictors in the model:
#'  \deqn{n_m = \frac{n}{1 - R^2}}{
#'        n[m] = n / (1 - R^2)}
#'
#' For a \bold{binomial} coefficient, the calculation uses
#' \eqn{P_0}{P[0]}, the probability given the null hypothesis and
#' \eqn{P_a}{P[a]}, the probability given the alternative hypothesis and
#' and the average probability \eqn{\bar{P} = \frac{P_0 + P_a}{2}}{
#'                                  Pbar = (P[0] + P[a]) /2}
#' The calculation is:
#'  \deqn{n = \frac{(z_{1-\alpha} \sqrt{2 \bar{P} (1 - \bar{P})} +
#'             z_{\code{beta}} \sqrt{P_0(1 - P_0) + P_a(1 - P_a)})^2}{
#'             (P_a + P_0)^2}}{
#'        n = (z[1-alpha](2Pbar(1 - Pbar)^0.5) + z[beta](P[0](1 - P[0]) + P[1](1 - P[1]))^0.5)^2 /
#'            (P[1] - P[0])^2}
#' An alternative given by Whitemore uses \eqn{\hat{P} = P(x=0)}{Phat = P(x=0)}.
#' \cr
#' The lead term in the equation below is used to correct for
#' large values of \eqn{\hat{P}}{Phat}:
#'  \deqn{n = (1 + 2P_0) \frac{(z_{1-\alpha} \sqrt{\frac{1}{1-\hat{P}} + \frac{1}{\hat{P}}} +
#'                             z_{\code{beta}} \sqrt{\frac{1}{1-\hat{P}} +
#'                                            \frac{1}{\hat{P} \exp{\hat{\beta}}}})^2}{
#'                       (P_0 \hat{\beta})^2}}{
#'        n = (1 + 2P[0]) * (z[1-alpha]sqrt(1/Phat + 1/(1+Phat)) +
#'                           z[beta]sqrt(1/Phat + 1/(Phat exp(Bhat))))^2 /
#'                           (P[0]Bhat)^2}
#' As above these can be adjusted in the multivariable case:
#'  \deqn{n_m = \frac{n}{1 - R^2}}{
#'        n[m] = n / (1 - R^2)}
#' In this case, Pearsons \eqn{R^2} correlation is between the
#' fitted values from a logistic regression with \code{coeff} as the response
#' and the other predictors as co-variates.
#' \cr
#' The calculation uses \eqn{\bar{P}}{Pbar}, the mean probability (mean of the
#' fitted values from the model):
#'  \deqn{R^2 = \frac{(\sum{i=1}^n (y_i - \bar{P})(P_i - \bar{P}))^2}{
#'                     \sum{i=1}^n (y_i - \bar{P})^2 \sum{i=1}^n (P_i - \bar{P})^2}}{
#'        R^2 = (sum(y[i] - Pbar)(P[i] - Pbar))^2 / 
#'              (sum(y[i] - Pbar)^2 * sum (P[i] - Pbar)^2)}
#' 
#' @note
#' The returned \code{list} has the additional
#' \code{class} of \code{"ss.glm"}.
#' \cr
#' The \code{print} method for this \code{class} does not 
#' show the attributes. 
#' 
#' @references Whitemore AS (1981).
#' Sample Size for Logistic Regression with Small Response Probability.
#' \emph{Journal of the American Statistical Association}. \bold{76}(373):27-32.
#' \href{http://dx.doi.org/10.2307/2287036}{
#'       JASA (paywall)}
#' \cr
#' JSTOR (free)
#' \cr
#' http://www.jstor.org/stable/2287036
#' 
#' @references Hsieh FY (1989).
#' Sample size tables for logistic regression.
#' \emph{Statistics in Medicine}. \bold{8}(7):795-802.
#' \href{http://dx.doi.org/10.1002/sim.4780080704}{
#'       Wiley (paywall)}.
#' \href{http://www.statpower.net/Content/312/Handout/Hsieh\%281989\%29.pdf}{
#'       statpower (free)}.
#' 
#' @references Fleiss J (2003).
#' \emph{Statistical methods for rates and proportions. 3rd ed}.
#' John Wiley, New York.
#' \href{dx.doi.org/10.1002/0471445428}{
#'       Wiley (paywall)}.
#' \href{http://books.google.com/books?isbn=1118625617}{
#'       Google books (free preview)}.
#'
#' @references Peduzzi P, Concato J, Kemper E,
#' Holford T R, Feinstein A R (1996).
#' A simulation study of the number of events
#' per variable in logistic regression analysis.
#' \emph{Journal of clinical epidemiology}. \bold{49}(12):1373-79.
#' \href{http://dx.doi.org/10.1016/S0895-4356(96)00236-3}{
#'       JCE (paywall)}.
#' \href{http://www.researchgate.net/publication/14236358}{
#'       ResearchGate (free)}.
#' 
#' @keywords htest
#' 
#' @examples
#' ## H&L 2nd ed. Section 8.5.
#' ## Results here are slightly different from the text due to rounding.
#' data(uis)
#' with(uis, prop.table(table(DFREE, TREAT), 2))
#' (g1 <- glm(DFREE ~ TREAT, data=uis, family=binomial))
#' ss(g1, coeff="TREATlong")
#' ## Pages 340 - 341.
#' ss(g1, coeff="TREATlong", OR=1.5, Px0=0.5)
#' ## standardize
#' uis <- within(uis, {
#'     AGES <- (AGE - 32) / 6
#'     NDRGTXS <- (NDRGTX - 5) / 5
#' })
#' ## Page 343.
#' ## results slightly different due to rounding
#' g1 <- glm(DFREE ~ AGES, data=uis, family=binomial) 
#' ss(g1, coeff="AGES", std=FALSE, OR=1.5)
#' ## Table 8.37. Page 344.
#' summary(g1 <- glm(DFREE ~ AGES + NDRGTXS + IVHX + RACE + TREAT,
#'                   data=uis, family=binomial))
#' ## Page 345.
#' ## results slightly different due to rounding
#' ss(g1, coeff="AGES", std=FALSE, OR=1.5)
#' ss(g1, coeff="TREATlong", std=FALSE, OR=1.5)
ss.glm <- function(x,
                   ..., 
                   alpha=0.05,
                   beta=0.8,
                   coeff=names(stats::coef(x))[2],
                   std=FALSE,
                   alternative=c("one.sided", "two.sided"),
                   OR=NULL,
                   Px0=NULL
                   ) {
    stopifnot(inherits(x, "glm"))
    stopifnot(x$family$link=="logit")
    stopifnot(!stats::is.empty.model(x))
    stopifnot(!is.na(coeff))
    stopifnot(coeff %in% colnames(stats::model.matrix(x)))
    x1 <- stats::model.matrix(x)[, coeff]
    if(length(unique(x1))==1) stop(
                 "Must be more than one value for the coefficient")
    alternative <- match.arg(alternative)
    if (alternative=="one.sided"){
        z1 <- stats::qnorm(1 - alpha)
    } else {
        z1 <- stats::qnorm(1 - alpha/2)
    }
    z2 <- stats::qnorm(beta)
    res1 <- vector(mode="list")
    ## check if multiple predictors
    multi1 <- length(attr(x$terms, "term.labels")) > 1
    ## other predictors in multivariable model
    if(multi1){
        pred1 <- stats::model.matrix(x)[,
                                        !colnames(
                                            stats::model.matrix(x))==coeff,
                                        drop=FALSE]
    }
### coefficient is continuous
    if (length(unique(x1)) > 2){
        if (std) x1 <- (x1 - mean(x1)) / stats::sd(x1)
        ## value of intercept term
        B0 <- stats::coef(x)[names(stats::coef(x))=="(Intercept)"]
        ## null probability
        P0 <- exp(B0) / (1 + exp(B0))
        if (is.null(OR)){
            B1 <- x$coefficients[names(x$coefficients)==coeff]
        } else {
            B1 <- log(OR)
        }
        P1 <- exp(B1) / (1 + exp(B1))
        ## Whitmore's method
        delta1 <- (1 + ((1 + B1^2) * exp((1.25 * B1)^2))) /
            (1 + exp((-0.25 * B1)^2))
        n1 <- (1 + 2 * P0 * delta1) * (z1 + (z2 * exp((-0.25 * B1^2))))^2 /
            (P0 * B1^2)
        if (multi1){
            ## get R^2, the multiple correlation between
            ## coeff and the other predictors in model
            R2 <- stats::summary.lm(stats::lm(x1 ~ 0 + pred1))$r.squared
            ## Hsieh's multivariable adjustment
            n1a <- n1 / (1 - R2)
            res1$ss <- c(n1, n1a)
            names(res1$ss) <- c("uni", "multi")
        } else {
            res1$ss <- n1
            names(res1$ss) <- c("uni")
        }
        attr(res1$ss, "interpret") <- c(
            "Sample size required to show",
            paste0("probability for ", coeff, " = ", format(P1)),
            paste0("rather than ", coeff," = ", format(P0)))
        attr(res1$ss, "alpha") <- alpha
        attr(res1$ss, "beta") <- beta
        attr(res1$ss, "alternative") <- alternative
    }
### binomial coefficient:
 if (length(unique(x1)) <= 2) {
     ## probability P(y=1|x=0)
     P0 <- sum(x$y[x1==0]) / length(x$y[x1==0])
     if (is.null(OR)){
         ## P(y=1|x=1)
         P1 <- sum(x$y[x1==1]) / length(x$y[x1==1])
     } else {
         P1 <- OR * P0 / ((1 - P0) + OR * P0)
     }
     Pbar <- (P0 + P1) / 2
     n1 <- (z1 * sqrt(2 * Pbar * (1 - Pbar)) +
            z2 * sqrt(P0 * (1 - P0) + P1 * (1 - P1)))^2 / ((P1 - P0)^2)
     ## alternative method, using sampling distribution of Wald statistic
     ## for the estimate of coefficient
     if (is.null(Px0)) Px0 <- sum(x1==0) / length(x$y)
     if (is.null(OR)){
         ## coefficient for predictor
         B1 <- x$coefficients[[2]]
     } else {
         B1 <- log(OR)
     }
     n2 <- (1 + 2 * P0) * (
         z1 * sqrt((1 / (1 - Px0)) + (1 / Px0)) +
         z2 * sqrt((1 / (1 - Px0)) + (1 / (Px0 * exp(B1)))))^2 / 
             (P0 * B1^2)
     if (multi1){
         g2 <- stats::glm(x1 ~ 0 + pred1, family=stats::binomial)
         ybar <- pbar <- mean(g2$fitted.values)
         R2 <- sum((g2$y - ybar) * (g2$fitted.values - pbar))^2 /
             ((sum((g2$y - ybar)^2)) * (sum((g2$fitted.values - pbar)^2)))
         ## Hsieh's multivariable adjustment
         n1a <- n1 / (1 - R2)
         n2a <- n2 / (1 - R2)
         res1$ss <- c(n1, round(n2/2), n1a, round(n2a/2))
         names(res1$ss) <- c("uni", "uni_alt", "multi", "multi_alt")
        } else {
            res1$ss <- c(n1, round(n2/2))
            names(res1$ss) <- c("uni", "uni_alt")
        }
      attr(res1$ss, "interpret") <- c(
         "Sample size required to show",
         paste0("P(y=1 | x=0) = ", format(P0)),
         paste0("is different from  P(y=1|x=1) = ", format(P1)))
 }
### events per covariate
    res1$epc <- min(sum(x$y), length(x$y) - sum(x$y) ) / (length(x$coefficients)-1)
    attr(res1$epc, "interpret") <- "Events per covariate; ideally >10"
###
    class(res1) <- c("ss.glm", class(res1))
    return(res1)
}
