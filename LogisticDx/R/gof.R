#' @name gof
#' @title Goodness of fit tests for binomial regression
#'
#' @rdname gof
#' @export
#'
gof <- function(x, ...){
    UseMethod("gof")
}
#' @rdname gof
#' @aliases gof.glm
#' @method gof glm
#' @export
#'
#' @include dx.R
#' @include genBinom.R
#'
#' @param x A regression model with class \code{glm} and
#' \code{x$family$family == "binomial"}.
#' @param ... Additional arguments when plotting the
#' receiver-operating curve. See:
#' \cr
#' ?pROC::roc
#' \cr
#' and
#' \cr
#' ?pROC::plot.roc
#' @param g Number of groups (quantiles) into which to
#' split observations for
#' the Hosmer-Lemeshow and the modified Hosmer-Lemeshow tests.
#' @param plotROC Plot a receiver operating curve?
#' 
#' @return A \code{list} of \code{data.table}s as follows:
#'
#' \item{ct}{Contingency table.}
#'
#' \item{chiSq}{\eqn{\chi^2}{Chi-squared} tests of the significance of the model.
#'   The tests are:
#'   \tabular{cl}{
#'     PrI \tab test of the Pearsons residuals, calculated by \emph{individual} \cr
#'     drI \tab test of the deviance residuals, calculated by \emph{individual} \cr
#'     PrG \tab test of the Pearsons residuals, calculated by covariate \emph{group} \cr
#'     drG \tab test of the deviance residuals, calculated by covariate \emph{group} \cr
#'     PrCT \tab test of the Pearsons residuals, calculated from the \emph{contingency table} \cr
#'     drCT \tab test of the deviance residuals, calculated from the \emph{contingency table} 
#' }}
#'
#' \item{ctHL}{\bold{C}ontingency \bold{t}able for the \bold{H}osmer-\bold{L}emeshow test.}
#'
#' \item{gof}{
#' Goodness-of-fit tests. These are:
#'  \itemize{
#'    \item HL Hosmer-Lemeshow's \eqn{C} statistic.
#'    \item mHL The modified Hosmer-Lemeshow test.
#'    \item OsRo Osius and Rojek's test of the link function.
#'    \item S Stukel's tests:
#'     \tabular{ll}{
#'      SstPgeq0.5 \tab score test for addition of vector \eqn{z1} \cr
#'      SstPl0.5 \tab score test for addition of vector \eqn{z2} \cr
#'      SstBoth  \tab score test for addition of vector \eqn{z2} \cr
#'      SllPgeq0.5  \tab log-likelihood test for
#'                       addition of vector \eqn{z1} \cr
#'      SllPl0.5  \tab log-likelihood test for
#'                     addition of vector \eqn{z2} \cr
#'      SllBoth \tab log-likelihood test for addition
#'             of vectors \eqn{z1} and \eqn{z2}
#'       }}}
#'
#' \item{R2}{R-squared like tests:
#'  \tabular{cl}{
#'   ssI \tab sum-of-squares, by \emph{individual} \cr
#'   ssG \tab sum-of-squares, by covariate \emph{group} \cr
#'   llI \tab log-likelihood, by \emph{individual} \cr
#'   llG \tab log-likelihood, by covariate \emph{group}.
#'  }}
#'
#' \item{auc}{
#'  Area under the receiver-operating curve (ROC)
#'  with 95 \% CIs.
#' }
#' Additionally, if \code{plotROC=TRUE}, a plot of the ROC.
#'
#' @details
#' Details of the elements in the returned \code{list} follow below:
#' \cr \cr
#' \bold{ct}:
#'  \cr
#'  A contingency table, similar to the output of \code{\link{dx}}.
#'  \cr
#'  The following are given per \emph{covariate group}:
#'   \tabular{cl}{
#'    n \tab number of observations \cr
#'    y1hat \tab predicted number of observations with \eqn{y=1} \cr
#'    y1 \tab actual number of observations with \eqn{y=1} \cr
#'    y0hat \tab predicted number of observations with \eqn{y=0} \cr
#'    y0 \tab actual number of observations with \eqn{y=0} 
#'   }
#' 
#' \bold{chiSq}:
#'  \cr
#'  \eqn{P \chi^2}{Chi-squared} tests of the significance
#'  of the model.
#'  \cr
#'  Pearsons test and the deviance \eqn{D} test are given.
#'  \cr
#'  These are calculated by indididual \code{I}, by covariate group \code{G}
#'  and also from the contingency table \code{CT} above.
#'  They are calculated as:
#'   \deqn{P \chi^2 = \sum_{i=1}^n Pr_i^2}{
#'         Chisq = SUM Pr[i]^2}
#'  and
#'   \deqn{D = \sum_{i=1}^n dr_i^2}{
#'         D = SUM dr[i]^2}
#'  The statistics should follow a
#'  \eqn{\chi^2}{chiSq} distribution with
#'  \eqn{n - p} degrees of freedom.
#'  \cr
#'  Here, \eqn{n} is the number of observations
#'  (taken individually or by covariate group) and \eqn{p}
#'  is the number
#'  pf predictors in the model.
#'  \cr
#'  A \bold{high} \eqn{p} value for the test suggests
#'  that the model is a poor fit.
#'  \cr
#'  The assumption of a \eqn{\chi^2}{chiSq} distribution
#'  is most valid when
#'  observations are considered by group.
#'  \cr
#'  The statistics from the contingency table should be
#'  similar to those obtained
#'  when caluclated by group.
#'  \cr
#'
#' \bold{ctHL}:
#'  \cr 
#'  The contingency table for the Hosmer-Lemeshow test.
#'  \cr
#'  The observations are ordered by probability, then
#'  grouped into \code{g} groups of approximately equal size.
#'  \cr
#'  The columns are:
#'   \tabular{cl}{
#'    P \tab the probability \cr
#'    y1 \tab the actual number of observations with \eqn{y=1} \cr
#'    y1hat \tab the predicted number of observations with \eqn{y=1} \cr
#'    y0 \tab the actual number of observations with \eqn{y=0} \cr
#'    y0hat \tab the predicted number of observations with \eqn{y=0} \cr
#'    n \tab the number of observations \cr
#'    Pbar \tab the mean probability, which is \eqn{\frac{nP}{\sum{n}}}{(n * P) / SUM(n)} 
#'   }
#'
#' \bold{gof}:
#'  \cr 
#'  All of these tests rely on assessing the effect of
#'  adding an additional variable to the model.
#'  \cr
#'  Thus a \bold{low} \eqn{p} value for any of these tests
#'  implies that the model is a \bold{poor} fit.
#' 
#'  \subsection{Hosmer and Lemeshow tests}{
#'   Hosmer and Lemeshows \eqn{C} statistic is based on:
#'   \eqn{y_k}{y[k]}, the number of observations
#'   where \eqn{y=1},
#'   \eqn{n_k}{n[k]}, the number of observations and
#'   \eqn{\bar{P_k}}{Pbar[k]}, the average probability
#'   in group \eqn{k}:
#'    \deqn{\bar{P_k} = \sum_{i=1}^{i=n_k} \frac{n_iP_i}{n_k}, \quad k=1,2, \ldots ,g}{
#'          Pbar[k] = SUM(i) n[i]P[i] / n[k], k=1,2...g}
#'   The statistic is:
#'    \deqn{C = \sum_{k=1}^g \frac{(y_k - n_k \bar{P_k})^2}{
#'                                 n_k \bar{P_k} (1 - \bar{P_k})}}{
#'          C = SUM(k=1...g) (y[k] - n[k]Pbar[k])^2 / n[k]Pbar[k](1-Pbar[k])}
#'   This should follow a \eqn{\chi^2}{chiSq} distribution
#'   with \code{g - 2} degrees of freedom.
#'   \cr \cr
#'   The \bold{modified} Hosmer and Lemeshow test is assesses the
#'   change in model deviance \eqn{D} when \code{G} is added as a predictor.
#'   That is, a linear model is fit as:
#'    \deqn{dr_i \sim G, \quad dr_i \equiv \mathrm{deviance residual}}{
#'          dr[i] ~ G, dr[i] = deviance residual}
#'   and the effect of adding \eqn{G} assessed with \code{anova(lm(dr ~ G))}. 
#'  }
#' 
#'  \subsection{Osius and Rojek's tests}{
#'   These are based on a \emph{power-divergence} statistic \eqn{PD_{\lambda}}{PD[l]} 
#'   (\eqn{\lambda = 1}{l=1} for Pearsons test) and 
#'   the standard deviation (herein, of a binomial distribution)
#'   \eqn{\sigma}{SD}.
#'   The statistic is:
#'    \deqn{Z_{OR} = \frac{PD_{\lambda} - \mu_{\lambda}}{\sigma_{\lambda}}}{
#'          Z[OR] = PD[l] - lbar / SD[l]}
#'   \cr
#'   For logistic regression, it is calculated as:
#'    \deqn{Z_{OR} = \frac{P \chi^2 - (n - p)}{\sqrt{2(n - \sum_{i=1}^n \frac{1}{n_i}) + RSS}}}{
#'          Z[OR] = (chiSq - (n - p)) / (2 * n * SUM 1/n[i])^0.5}
#'   where \eqn{RSS} is the residual sum-of-squares from
#'   a weighted linear regression:
#'    \deqn{ \frac{1 - 2P_i}{\sigma_i} \sim X, \quad \mathrm{weights=} \sigma_i}{
#'          (1 - 2 * P[i]) / SD[i] ~ X, weights = SD[i]}
#'   Here \eqn{\bold{X}} is the matrix of model predictors.
#'   \cr
#'   A two-tailed test against a standard normal distribution \eqn{N ~ (0, 1)}
#'   should \emph{not} be significant.
#'   }
#' 
#'  \subsection{Stukels tests}{
#'   These are based on the addition of the vectors:
#'    \deqn{z_1 = \mathrm{Pgeq0.5} = \mathrm{sign}(P_i \geq 0.5)}{
#'          z[1] = Pgeq0.5 = sign(P[i] >= 0.5)}
#'   and / or
#'    \deqn{z_2 = \mathrm{Pl0.5} = \mathrm{sign}(P_i < 0.5)}{
#'          z[2] = Pl0.5 = sign(P[i] < 0.5)}
#'  to the existing model predictors.
#'  \cr
#'  The model fit is compared to the original using the
#'  score (e.g. \code{SstPgeq0.5}) and likelihood-ratio (e.g. \code{SllPl0.5}) tests.
#'  These models should \emph{not} be a significantly better fit to the data.
#'  }
#'
#' \bold{R2}:
#'  \cr \cr
#'  Pseudo-\eqn{R^2} comparisons of the predicted values from the
#'  fitted model vs. an intercept-only model.
#' 
#'  \subsection{sum-of-squares}{
#'   The sum-of-squres (linear-regression) measure based on
#'   the squared Pearson correlation coefficient
#'   by \emph{individual} is based on the mean probability:
#'    \deqn{\bar{P} = \frac{\sum n_i}{n}}{
#'         Pbar = SUM(n[i]) / n}
#'   and is given by:
#'    \deqn{R^2_{ssI} = 1 - \frac{\sum (y_i - P_i)^2}{\sum (y_i - \bar{P})^2}}{
#'         R2[ssI] = 1 - SUM(y[i] - P[i])^2 / SUM(y[i] - Pbar)^2}
#'   The same measure, by \emph{covariate group}, is:
#'    \deqn{R^2_{ssG} = 1 - \frac{\sum (y_i - n_iP_i)^2}{\sum (y_i - n_i\bar{P})^2}}{
#'          R2[ssG] = 1 - SUM(y[i] - n[i] * P[i])^2 / SUM(y[i] - n[i] * Pbar)^2}
#'  }
#'
#'  \subsection{log-likelihood}{
#'   The log-likelihood based \eqn{R^2} measure per \emph{individual} is
#'   based on:
#'    \itemize{
#'     \item \eqn{ll_0}{ll[0]}, the log-likelihood of the intercept-only model
#'     \item \eqn{ll_p}{ll[p]}, the log-likelihood of the model with \eqn{p} covariates
#'   }
#'   It is calculated as 
#'    \deqn{R^2_{llI} = \frac{ll_0 - ll_p}{ll_0} = 1 - \frac{ll_p}{ll_0}}{
#'          R2[llI] = (ll[0] - ll[p]) / ll[0]}
#'   This measure per \emph{covariate group} is based on
#'   \eqn{ll_s}{ll[s]}, the log-likelihood for the \emph{saturated} model,
#'   which is calculated from the model deviance \eqn{D}:
#'    \deqn{ll_s = ll_p - \frac{D}{2}}{
#'          ll[s] = ll[p] - D / 2}
#'   It is cacluated as:
#'    \deqn{R^2_{llG} = \frac{ll_0 - ll_p}{ll_0 - ll_s}}{
#'          R2[llG] = (ll[0] - ll[p]) / (ll[0] - ll[s])}
#'  }
#' 
#' \bold{auc}:
#'  \cr \cr
#'  The area under the receiver-operating curve.
#'  \cr
#'  This may broadly be interpreted as:
#'   \tabular{cc}{
#'    auc \tab Discrimination \cr
#'    \eqn{\mathrm{auc} = 0.5}{auc=0.5} \tab useless \cr
#'    \eqn{0.7 \leq \mathrm{auc} < 0.8}{0.7 <= auc < 0.8}
#'      \tab acceptable \cr
#'    \eqn{0.8 \leq \mathrm{auc} < 0.9}{0.8 <= auc < 0.9}
#'      \tab excellent \cr
#'    \eqn{\mathrm{auc} \geq 0.9}{auc >= 0.9}
#'      \tab outstanding
#'   }
#'  \eqn{\mathrm{auc} \geq 0.9}{auc >= 0.9} occurs rarely
#'  as this reuqires almost
#'  complete separation/ perfect classification.
#'  \cr
#'
#'
#' @note
#' The returned \code{list} has the additional
#' \code{class} of \code{"gof.glm"}.
#' \cr
#' The \code{print} method for this \code{class} shows
#' \emph{only} those results
#' which have a \eqn{p} value.
#'
#' @references Osius G & Rojek D, 1992.
#' Normal goodness-of-fit tests for multinomial models
#' with large degrees of freedom.
#' \emph{Journal of the American Statistical Association}.
#' \bold{87}(420):1145-52.
#' \href{http://dx.doi.org/10.1080/01621459.1992.10476271}{
#'       JASA (paywall)}.
#' \cr
#' JSTOR (free):
#' \cr
#' http://www.jstor.org/stable/2290653
#'
#' @references Hosmer D, Hosmer T, Le Cessie S & Lemeshow S (1997).
#' A comparison of goodness-of-fit tests for
#' the logistic regression model.
#' \emph{Statistics in Medicine}. \bold{16}(9):965-80.
#' \href{dx.doi.org/10.1002/(SICI)1097-0258(19970515)16:9<965::AID-SIM509>3.0.CO;2-O}{
#'       Wiley (paywall)}.
#' \href{https://stat.duke.edu/~zo2/dropbox/goflogistic.pdf}{
#'       Duke University (free)}.
#'
#' @references Mittlboch M, Schemper M (1996).
#' Explained variation for logistic regression.
#' \emph{Statistics in Medicine}. \bold{15}(19):1987-97.
#' \href{dx.doi.org/10.1002/(SICI)1097-0258(19961015)15:19<1987::AID-SIM318>3.0.CO;2-9}{
#'       Wiley (paywall)}.
#' \href{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.477.3328&rep=rep1&type=pdf}{
#'       CiteSeerX / Penn State University (free)}.
#'
#' @author Modified Hosmer & Lemeshow goodness of fit test:
#' adapted from existing work by Yongmei Ni.
#' 
#' @keywords htest
#' 
#' @examples
#' ## H&L 2nd ed. Sections 5.2.2, 5.2.4, 5.2.5. Pages 147-167.
#' data(uis)
#' uis <- within(uis, {
#'     NDRGFP1 <- 10 / (NDRGTX + 1)
#'     NDRGFP2 <- NDRGFP1 * log((NDRGTX + 1) / 10)
#' })
#' g1 <- gof(glm(DFREE ~ AGE + NDRGFP1 + NDRGFP2 + IVHX +
#'               RACE + TREAT + SITE +
#'               AGE:NDRGFP1 + RACE:SITE,
#'               family=binomial, data=uis), plot=FALSE)
#' g1
#' unclass(g1)
#' attributes(g1$gof)
gof.glm <- function(x, ..., g=10, plotROC=TRUE){
    stopifnot(inherits(x, "glm"))
    stopifnot(x$family$family=="binomial")
    ## for R CMD check
    Pr <- dr <- n <- P <- yhat <- y <- NULL
    y1 <- y1hat <- y0 <- y0hat <- NULL
    chiSq <- G <- Pbar <- df <- pVal <- Z <- NULL
    ### hold results
    res1 <- vector(mode="list")
    attr(res1, "link") <- x$family$link
###
### chi-square tests
###
    ## by observation
    ## Pearson residuals
    PrI <- sum(stats::residuals.glm(x, type="pearson")^2)
    ## deviance residuals
    ## same as Residual Deviance from summary(model)
    drI <- sum(stats::residuals.glm(x, type="deviance")^2)
    degfI <- stats::df.residual(x)
    ## by covariate group
    dx1 <- dx(x)
    ## Pearson residuals
    PrG <- dx1[, sum(Pr^2)]
    ## deviance residuals
    drG <- dx1[, sum(dr^2)]
    ## degrees of freedom =
    ## no. covariate patterns (with any observations)
    ##  - (no. coefficients)
    degfG <- nrow(dx1) - length(x$coefficients)
### the above doesn't work well as nrow(dx1) approaches nrow(model$data)
### so also use a contingency table
    dx2 <- dx1[, list(n, P, yhat, y, n * (1 - P), n - y)]
    data.table::setnames(dx2, c("n", "P", "y1hat", "y1", "y0hat", "y0"))
### manual chi-sq test
    PrCT <- dx2[, sum(((y1 - y1hat)^2) / y1hat) + sum(((y0 - y0hat)^2) / y0hat)]
    drCT <- dx2[, 2 * sum(y1 * log(y1 / y1hat), na.rm=TRUE)] +
        dx2[, 2 * sum(y0 * log(y0 / y0hat), na.rm=TRUE)]
    degfCT <- nrow(dx2) - length(x$coefficients)
###
    res1$ct <- dx2[, list(n, y1hat, y1, y0hat, y0)]
    res1$chiSq <- data.table::data.table(
        "test"=c("PrI", "drI", "PrG", "drG", "PrCT", "drCT"),
        "chiSq"=c(PrI, drI, PrG, drG, PrCT, drCT),
        "df"=rep(c(degfI, degfG, degfCT), each=2L))
    res1$chiSq[, "pVal" := stats::pchisq(chiSq, df, lower.tail=FALSE)]
    attr(res1$chiSq, "interpret") <- c(
        "Low p value: reject H0 that coefficients are poor predictors.",
        "I.e. coefficients are significant predictors.")
###
### GOF tests
###
### Hosmer-Lemeshow
    ## sort by probability
    dx1 <- dx1[order(P), ]
    ## base group size
    gs1 <- dx1[, sum(n) %/% g]
    g1 <- rep(gs1, g)
    ## remainer
    mod1 <- dx1[, sum(n) %% g]
    g1[seq(1, g-1, by=g/mod1)] <- gs1 + 1
    dx1[, "G" := cut(cumsum(n), breaks=c(0, cumsum(g1)))]
    dx1[, "G" := factor(G,
                  labels=dx1[, format(max(P), digits=3), by=G]$V1)]
    dx2 <- dx2[order(P), ]
    dx2[, "G" := dx1[, G]]
    res1$ctHL <- dx2[, list(sum(y1), sum(y1hat), sum(y0),
                            sum(y0hat), sum(n), sum(n*P/length(n))), by=G]
    res1$ctHL <- dx2[, list(sum(y1), sum(y1hat), sum(y0),
                            sum(y0hat), sum(n), sum(n * P / sum(n))), by=G]
    data.table::setnames(res1$ctHL, c("P", "y1", "y1hat", "y0", "y0hat", "n", "Pbar"))
    ## average probability by group
    ## dx1[, "Pg" := sum(n * P / sum(n)), by=G]
    HL1 <- sum(res1$ctHL[, (y1 - n * Pbar)^2 / (n * Pbar * (1 - Pbar))])
    ## modified H&L
    HL2 <- dx1[, stats::anova(stats::lm(dr ~ G))][, c("Df", "F value", "Pr(>F)")]["G", ]
### Osius & Rojek test
    ## A1 = correction factor for the variance
    A1 <- dx1[, 2 * (nrow(dx1) - sum(1 / n))]
    ## weighted least-squares regression
    l1 <- length(x$coefficients)
    lm1 <- stats::lm(dx1[, (1 - 2 * P) / (n * P * (1 - P))] ~
              as.matrix(dx1[, seq.int(l1), with=FALSE]),
              weights=dx1[, (n * P * (1 - P))])
    ## residual sum of squares
    RSS1 <- sum(stats::residuals(lm1, type="pearson")^2)
    ## PdJ = chi-squared test for model, by covariate group
    z1 <- (PrG - (nrow(dx1) - l1 )) / sqrt(A1 + RSS1)
### Stukels test
    link1 <- x$family$linkfun(x$fitted.values)
    ## probability >= 0.5
    Pgre0.5 <- 0.5 * link1^2 * (x$fitted.values >= 0.5)
    ## probability < 0.5
    Pl0.5 <- -0.5 * link1^2 * (x$fitted.values < 0.5)
    Z1 <- abs(statmod::glm.scoretest(x, matrix(c(Pgre0.5, Pl0.5), ncol=2)))
    chi1 <- sum(Z1^2)
    ##
    environment(x$formula) <-  environment()
    m1 <- stats::update(x, formula=stats::update.formula(x$formula, ~ . + Pgre0.5))
    LR1 <- -2 * (x$deviance / -2 - m1$deviance / -2)
    m1 <- stats::update(x, formula=stats::update.formula(x$formula, ~ . + Pl0.5))
    LR2 <- -2 * (x$deviance / -2 - m1$deviance / -2)
    m1 <- stats::update(x, formula=stats::update.formula(x$formula, ~ . + Pgre0.5 + Pl0.5))
    LR3 <- -2 * (x$deviance / -2 - m1$deviance / -2)
    ##
    res1$gof <- data.table::data.table(
        "test"=c("HL","mHL", "OsRo",
            "SstPgeq0.5", "SstPl0.5", "SstBoth",
            "SllPgeq0.5", "SllPl0.5", "SllBoth"),
        "stat"=c("chiSq", "F", "Z", "Z", "Z", "chiSq", "chiSq", "chiSq", "chiSq"),
        "val"=c(HL1, HL2[[2]], z1, Z1, chi1, LR1, LR2, LR3),
        "df"=c(g-2, HL2[[1]], NA, NA, NA, 2, 1, 1, 2),
        "pVal"=c(stats::pchisq(HL1, g - 2, lower.tail=FALSE),
            HL2[[3]],
            stats::pnorm(abs(z1), lower.tail=FALSE) * 2,
            stats::pnorm(abs(Z1[1]), lower.tail=FALSE) * 2,
            stats::pnorm(abs(Z1[2]), lower.tail=FALSE) * 2,
            stats::pchisq(chi1, 2, lower.tail=FALSE),
            stats::pchisq(LR1, 1, lower.tail=FALSE),
            stats::pchisq(LR2, 1, lower.tail=FALSE),
            stats::pchisq(LR3, 2, lower.tail=FALSE)))
    attr(res1$gof, "interpret") <- c(
        "Low p value: reject H0 that the model is a good fit.",
        "I.e. model is a poor fit.")
    attr(res1$gof, "g") <- paste0(g, " groups")
###
### pseudo R^2 tests
###
    ## linear regression-like sum of squares R^2
    ## Pearson r^2 correlation of observed outcome with predicted probability
    ssI <- 1 - (sum((x$y - x$fitted.values)^2)) / (sum((x$y - mean(x$fitted.values))^2))
    ## using covariate patterns
    ssG <- 1 - dx1[, sum((y - n * P)^2) / sum((y- n * mean(P))^2)]
    ## log-likelihood-based R^2
    ## stopifnot(logLik(x)==(x$deviance / -2))
    llI <- 1 - ((x$deviance / -2) / (x$null.deviance / -2))
    ## using covariate patterns
    ## deviance residual from covariate table above
    llG <- ((x$null.deviance / -2) - (x$deviance / -2)) /
        ((x$null.deviance / -2) - ((x$deviance / -2) + (0.5 * drG)))
    res1$R2 <- data.table::data.table(
        "method"=c("ssI", "ssG", "llI", "llG"),
        "R2"=c(ssI, ssG, llI, llG))
    attr(res1$R2, "details") <- c(
        "ss = linear regression-like sum of squares R^2",
        "ll = log-likelihood based R^2")
###
### ROC curve
    ## plot
    ## get name of y/ outcome variable from formula
    r1 <- pROC::roc(response=x$y,
                    predictor=x$fitted.values,
                    ci=TRUE, percent=TRUE, ...)
    ## 0.5 = chance, aim >0.7
    if(plotROC){
        pROC::plot.roc(r1, print.auc=TRUE, grid=TRUE,
                       print.auc.cex=0.8,
                       main="Receiver Operating Curve")
    }
    ## AUC and CI
    res1$auc <- c(as.vector(r1$auc),
                  as.vector(r1$ci)[c(1, 3)])
    names(res1$auc) <- c("auc",
                         paste0("lower ",
                                100 * attr(r1$ci, "conf.level"),
                                "% CI"),
                         paste0("upper ",
                                100 * attr(r1$ci, "conf.level"),
                                "% CI"))
    attr(res1$auc, "interpret") <- c(
        "auc = 0.5       --> useless",
        "0.7 < auc < 0.8 --> good",
        "0.8 < auc < 0.9 --> excellent")
###
    class(res1) <- c("gof.glm", class(res1))
    return(res1)
}
