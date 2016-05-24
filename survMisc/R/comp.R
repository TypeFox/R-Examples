#' @name comp
#' @title compare survival curves
#' 
#' @include ten.R
#' @include asWide.R
#' @include COV.R
#' @include predict.R
#' @include sf.R
#' @include print.R
#' 
#' @rdname comp
#' @export comp
#'
comp <- function(x, ...) UseMethod("comp")
#'
#' @param x A \code{tne} object
#' @param p \eqn{p} for Fleming-Harrington test
#' @param q \eqn{q} for Fleming-Harrington test
#' @param scores scores for tests for trend
#' @inheritParams sf.ten
#' 
#' @return The \code{tne} object is given
#'  additional \code{attributes}.
#'  \cr
#' The following are always added:
#' \item{lrt}{The \bold{l}og-\bold{r}ank family of \bold{t}ests}
#' \item{lrw}{The \bold{l}og-\bold{r}ank \bold{w}eights (used in calculating the tests).}
#' An additional item depends on the number of covariate groups.
#' \cr
#' If this is \eqn{=2}:
#' \item{sup}{The \bold{sup}remum or Renyi family of tests}
#' and if this is \eqn{>2}:
#' \item{tft}{Tests for trend. This is given as a \code{list}, 
#'  with the statistics and the scores used.}
#' 
#' @details
#' The \bold{log-rank} tests are formed from the following elements,
#' with values for each time where there is at least one event:
#' \itemize{
#'  \item \eqn{W_i}{W[i]}, the weights, given below.
#'  \item \eqn{e_i}{e[i]}, the number of events (per time).
#'  \item \eqn{\hat{e_i}}{P[i]}, the number of \emph{predicted} events,
#'                               given by \code{\link{predict}}.
#'  \item \eqn{COV_i}{COV[, , i]}, the covariance matrix for time \eqn{i},
#'                                 given by \code{\link{COV}}.
#' }
#' It is calculated as:
#'  \deqn{Q_i = \sum{W_i (e_i - \hat{e}_i)}^T
#'              \sum{W_i \hat{COV_i} W_i^{-1}}
#'              \sum{W_i (e_i - \hat{e}_i)}}{
#'       Q[i] = sum(W[i] * (e[i] - P[i]))^T *
#'              sum(W[i] * COV[, , i] * W[i])^-1 *
#'              sum(W[i] * (e[i] - P[i]))}
#' 
#' If there are \eqn{K} groups, then \eqn{K-1} are selected (arbitrary).
#'  \cr
#' Likewise the corresponding variance-covariance matrix is reduced to the
#' appropriate \eqn{K-1 \times K-1}{K-1 * K-1} dimensions.
#'  \cr
#' \eqn{Q} is distributed as chi-square with \eqn{K-1} degrees of freedom.
#' \cr \cr
#' For \eqn{2} covariate groups, we can use:
#' \itemize{
#'  \item \eqn{e_i}{e[i]} the number of events (per time).
#'  \item \eqn{n_i}{e[i]} the number at risk overall.
#'  \item \eqn{e1_i}{e1[i]} the number of events in group \eqn{1}.
#'  \item \eqn{n1_i}{n1[i]} the number at risk in group \eqn{1}.
#' } 
#' Then:
#'  \deqn{Q = \frac{\sum{W_i [e1_i - n1_i (\frac{e_i}{n_i})]} }{
#'                 \sqrt{\sum{W_i^2 \frac{n1_i}{n_i}
#'                            (1 - \frac{n1_i}{n_i})
#'                            (\frac{n_i - e_i}{n_i - 1}) e_i }}}}{
#'        Q = sum(W[i] * (e1[i] - n1[i] * e[i] / n[i])) /
#'           sqrt(sum(W[i]^2 * e1[i] / e[i] * (1 - n1[i] / n[i]) * (n[i] - e[i] / (n[i] - 1)) *e[i]))}
#' Below, for the Fleming-Harrington weights, 
#' \eqn{\hat{S}(t)}{S(t)} is the Kaplan-Meier (product-limit) estimator.
#'  \cr
#' Note that both \eqn{p} and \eqn{q} need to be \eqn{\geq 0}{>=0}.
#'  \cr \cr
#' The weights are given as follows:
#' \tabular{cll}{
#'  \eqn{1} \tab log-rank \tab \cr
#'  \eqn{n_i}{n[i]} \tab Gehan-Breslow generalized Wilcoxon \tab \cr
#'  \eqn{\sqrt{n_i}}{sqrt(n[i])} \tab Tarone-Ware \tab \cr
#'  \eqn{S1_i}{S1[i]} \tab Peto-Peto's modified survival estimate \tab
#'                    \eqn{\bar{S}(t)=\prod{1 - \frac{e_i}{n_i + 1}}}{
#'                             S1(t) = cumprod(1 - e / (n + 1))} \cr
#'  \eqn{S2_i}{S2[i]} \tab modified Peto-Peto (by Andersen) \tab
#'                     \eqn{\tilde{S}(t)=\bar{S} - \frac{n_i}{n_i + 1}}{
#'                               S2(t) = S1[i] * n[i] / (n[i] + 1) } \cr
#'  \eqn{FH_i}{FH[i]} \tab Fleming-Harrington \tab
#'                   The weight at \eqn{t_0 = 1} and thereafter is:
#'                    \eqn{\hat{S}(t_{i-1})^p [1-\hat{S}(t_{i-1})^q]}{
#'                         S(t[i - 1])^p * (1 - S(t)[i - 1]^q)}
#' }
#' The \bold{supremum (Renyi)} family of tests are designed
#'  to detect differences in survival curves which \emph{cross}.
#' \cr
#' That is, an early difference in survival in favor of one group
#'  is balanced by a later reversal.
#' \cr
#' The same weights as above are used.
#' \cr
#' They are calculated by finding
#' \deqn{Z(t_i) = \sum_{t_k \leq t_i} W(t_k)[e1_k - n1_k\frac{e_k}{n_k}], \quad i=1,2,...,k}{
#'       Z(t[i]) = SUM W(t[k]) [ e1[k] - n1[k]e[k]/n[k] ]}
#' (which is similar to the numerator used to find \eqn{Q}
#' in the log-rank test for 2 groups above).
#' \cr
#' and it's variance:
#' \deqn{\sigma^2(\tau) = \sum_{t_k \leq \tau} W(t_k)^2 \frac{n1_k n2_k (n_k-e_k) e_k}{n_k^2 (n_k-1)} }{
#'       simga^2(tau) = sum(k=1, 2, ..., tau) W(t[k]) (n1[k] * n2[k] * (n[k] - e[k]) * 
#'                                              e[k] / n[k]^2 * (n[k] - 1) ] }
#' where \eqn{\tau}{tau} is the largest \eqn{t}
#' where both groups have at least one subject at risk.
#' \cr \cr
#' Then calculate:
#' \deqn{ Q = \frac{ \sup{|Z(t)|}}{\sigma(\tau)}, \quad t<\tau }{
#' Q = sup( |Z(t)| ) / sigma(tau), t<tau}
#' When the null hypothesis is true,
#' the distribution of \eqn{Q} is approximately
#'  \deqn{Q \sim \sup{|B(x)|, \quad 0 \leq x \leq 1}}{
#'        Q ~ sup( |B(x)|, 0 <= x <= 1)}
#' And for a standard Brownian motion (Wiener) process:
#'  \deqn{Pr[\sup|B(t)|>x] = 1 - \frac{4}{\pi}
#'                           \sum_{k=0}^{\infty}
#'                           \frac{(- 1)^k}{2k + 1} \exp{\frac{-\pi^2(2k + 1)^2}{8x^2}}}{
#'        Pr[sup|B(t)|>x] = 1 - 4/pi sum((-1)^k / (2 * k + 1) * exp(-pi^2 (2k + 1)^2 / x^2))}
#' \bold{Tests for trend} are designed to detect ordered differences in survival curves.
#' \cr
#' That is, for at least one group:
#' \deqn{S_1(t) \geq S_2(t) \geq ... \geq S_K(t) \quad t \leq \tau}{
#'       S1(t) >= S2(t) >= ... >= SK(t) for t <= tau}
#' where \eqn{\tau}{tau} is the largest \eqn{t} where all groups have at least one subject at risk.
#' The null hypothesis is that
#' \deqn{S_1(t) = S_2(t) = ... = S_K(t) \quad t \leq \tau}{
#'       S1(t) = S2(t) = ... = SK(t) for t <= tau}
#' Scores used to construct the test are typically \eqn{s = 1,2,...,K},
#' but may be given as a vector representing a numeric characteristic of the group.
#' \cr
#' They are calculated by finding:
#' \deqn{ Z_j(t_i) = \sum_{t_i \leq \tau} W(t_i)[e_{ji} - n_{ji} \frac{e_i}{n_i}], 
#'                    \quad j=1,2,...,K}{
#'       Z[t(i)] = sum(W[t(i)] * (e[j](i) - n[j](i) * e(i) / n(i)))}
#' The test statistic is:
#' \deqn{Z = \frac{ \sum_{j=1}^K s_jZ_j(\tau)}{\sqrt{\sum_{j=1}^K \sum_{g=1}^K s_js_g \sigma_{jg}}} }{
#'       Z = sum(j=1, ..., K) s[j] * Z[j] / sum(j=1, ..., K) sum(g=1, ..., K) 
#'                                          s[j] * s[g] * sigma[jg]}
#' where \eqn{\sigma}{sigma} is the the appropriate element in the
#' variance-covariance matrix (see \code{\link{COV}}).
#' \cr
#' If ordering is present, the statistic \eqn{Z} will be greater than the 
#' upper \eqn{\alpha}{alpha}-th
#' percentile of a standard normal distribution.
#'
#' @note Regarding the Fleming-Harrington weights: 
#' \itemize{
#'  \item \eqn{p = q = 0} gives the log-rank test, i.e. \eqn{W=1}
#'  \item \eqn{p=1, q=0} gives a version of the Mann-Whitney-Wilcoxon test
#'   (tests if populations distributions are identical)
#'  \item \eqn{p=0, q>0} gives more weight to differences later on
#'  \item \eqn{p>0, q=0} gives more weight to differences early on
#' }
#' The example using \code{alloauto} data illustrates this.
#' Here the log-rank statistic
#' has a p-value of  around 0.5
#' as the late advantage of allogenic transplants
#' is offset by the high early mortality. However using
#' Fleming-Harrington weights of \eqn{p=0, q=0.5},
#' emphasising differences later in time, gives a p-value of 0.04.
#'  \cr
#' Stratified models (\code{stratTen}) are \emph{not} yet supported.
#'
#' @references Gehan A.
#' A Generalized Wilcoxon Test for Comparing Arbitrarily
#' Singly-Censored Samples.
#' Biometrika 1965 Jun. 52(1/2):203--23.
#' \href{http://www.jstor.org/stable/2333825}{JSTOR}
#' @references Tarone RE, Ware J 1977
#' On Distribution-Free Tests for Equality of Survival Distributions.
#' \emph{Biometrika};\bold{64}(1):156--60.
#' \href{http://www.jstor.org/stable/2335790}{JSTOR}
#' @references Peto R, Peto J 1972
#' Asymptotically Efficient Rank Invariant Test Procedures.
#' \emph{J Royal Statistical Society} \bold{135}(2):186--207.
#' \href{http://www.jstor.org/stable/2344317}{JSTOR}
#' @references Fleming TR, Harrington DP, O'Sullivan M 1987
#' Supremum Versions of the Log-Rank and Generalized Wilcoxon Statistics.
#' \emph{J  American Statistical Association} \bold{82}(397):312--20.
#' \href{http://www.jstor.org/stable/2289169}{JSTOR}
#' @references Billingsly P 1999
#' \emph{Convergence of Probability Measures.}
#' New York: John Wiley & Sons.
#' \href{http://dx.doi.org/10.1002/9780470316962}{Wiley (paywall)}
#' 
#' @examples
#' ## Two covariate groups
#' ## K&M 2nd ed. Example 7.2, Table 7.2, pp 209--210.
#' data("kidney", package="KMsurv")
#' t1 <- ten(Surv(time=time, event=delta) ~ type, data=kidney)
#' comp(t1, p=c(0, 1, 1, 0.5, 0.5), q=c(1, 0, 1, 0.5, 2))
#' ## see the weights used
#' attributes(t1)$lrw
#' ## supremum (Renyi) test; two-sided; two covariate groups
#' ## K&M 2nd ed. Example 7.9, pp 223--226.
#' data("gastric", package="survMisc")
#' g1 <- ten(Surv(time, event) ~ group, data=gastric)
#' comp(g1)
#' ## Three covariate groups
#' ## K&M 2nd ed. Example 7.4, pp 212-214.
#' data("bmt", package="KMsurv")
#' b1 <- ten(Surv(time=t2, event=d3) ~ group, data=bmt)
#' comp(b1, p=c(1, 0, 1), q=c(0, 1, 1))
#' ## Tests for trend
#' ## K&M 2nd ed. Example 7.6, pp 217-218.
#' data("larynx", package="KMsurv")
#' l1 <- ten(Surv(time, delta) ~ stage, data=larynx)
#' comp(l1)
#' attr(l1, "tft")
#' ### see effect of F-H test
#' data("alloauto", package="KMsurv")
#' a1 <- ten(Surv(time, delta) ~ type, data=alloauto)
#' comp(a1, p=c(0, 1), q=c(1, 1))
#' 
#' @rdname comp
#' @aliases comp.ten
#' @method comp ten
#' @export
#'
comp.ten <- function(x,
                     ...,
                     p=1,
                     q=1,
                     scores=seq.int(attr(x, "ncg")),
                     reCalc=FALSE){
    if (!reCalc & !is.null(attr(x, "lrt"))) {
        print(attr(x, "lrt"))
        print(if (!is.null(attr(x, "sup"))) attr(x, "sup") else attr(x, "tft"))
        return(invisible())
    }
    stopifnot(attr(x, "ncg") >= 2)
    stopifnot(length(p)==length(q))
    ## number of F-H tests
    fh1 <- length(p)
    if (!attr(x, "sorted")=="t") data.table::setkey(x, t)
    ## times with at least one event
    t1 <- x[e>0, t, by=t][, t]
    ## WEIGHTS
    wt1 <- data.table::data.table(
        array(data=1, dim=c(length(t1), 5L + fh1)))
    ## names for F-H tests
    FHn <- paste("FH_p=", p, "_q=", q, sep="")
    ## names for weights
    n1 <- c("1", "n", "sqrtN", "S1", "S2", FHn)
    data.table::setnames(wt1, n1)
    ## Gehan-Breslow generalized Wilcoxon, weight = n
    data.table::set(wt1, j="n",
                    value=x[e>0, max(n), by=t][, V1])
    ## Tarone-Ware, weight = sqrt(n)
    data.table::set(wt1, j="sqrtN",
                    value=wt1[, sqrt(.SD), .SDcols="n"])
    ## Peto-Peto, weight = S(t) = modified estimator of survival function
    data.table::set(wt1, j="S1",
                    value=cumprod(x[e > 0, 1 - sum(e) / (max(n) + 1), by=t][, V1]))
    ## modified Peto-Peto (by Andersen), weight = S(t)n / n+1
    data.table::set(wt1, j="S2",
                    value=wt1[, S1] * x[e > 0, max(n) / (max(n) + 1), by=t][, V1])
    ## Fleming-Harrington
    S3 <- sf(x=x[e>0, sum(e), by=t][, V1], n=x[e>0, max(n), by=t][, V1], what="S")
    ## weight of first 1st element is 1 as depends on [i-1]
    S3 <- c(1, S3[seq.int(length(S3) - 1L)])
    ## assign to names
    ##  SIMPLIFY = FALSE returns list rather than matrix
    wt1[, (FHn) := mapply(function(p, q) S3^p * ((1 - S3)^q),
                          p, q, SIMPLIFY=FALSE)]
    ## Hold results:
    ##  if 2 groups res1 = log-rank tests
    ##  if >2 groups, res1 = tests for trend
    n2 <- c("W", "Q", "Var", "Z",
            "pNorm", "chiSq", "df", "pChisq")
    res1 <- data.table::data.table(matrix(0, nrow=ncol(wt1), ncol=length(n2)))
    data.table::setnames(res1, n2)
    data.table::set(res1, j=1L, value=n1)
    predict(x)
    ## events minus predicted
    eMP1 <- attr(x, "pred")
    eMP1 <- eMP1[rowSums(eMP1) > 0, ]
    ## covariance
    COV(x)
    cov1 <- attr(x, "COV")
    if (is.null(dim(cov1))) {
        cov1 <- cov1[names(cov1) %in% t1]
    } else {
        ## 3rd dimension = times
        cov1 <- cov1[, , dimnames(cov1)[[3]] %in% t1]
    }
    ## number of covariate groups
    ncg1 <- attr(x, "ncg")
### 2 groups only:
    if (ncg1==2) {
### log-rank family
        ## make observed - expected for one group
        eMP1 <- unlist(eMP1[, .SD, .SDcols=length(eMP1)])
        data.table::set(res1, j="Q",
                        value=colSums(wt1 * eMP1))
        data.table::set(res1, j="Var",
                        value=colSums(wt1^2 * cov1))
### supremum tests
###  aka Renyi statistics
###  (analagous to 2-sample Kolmogorov-Smirnov test)
        n3 <- c("W", "maxAbsZ", "Var", "Q", "pSupBr")
        res2 <- data.table::data.table(matrix(0, nrow=5 + fh1, ncol=length(n3)))
        data.table::setnames(res2, n3)
        data.table::set(res2, j=1L, value=n1)
        data.table::set(res2, j="maxAbsZ",
                        value=sapply(abs(cumsum(eMP1 * wt1)), max))
        data.table::set(res2, j="Var",
                        value=res1[, Var])
        res2[, "Q" := maxAbsZ / sqrt(Var)]
        res2[, "pSupBr" := sapply(Q, probSupBr)]
        data.table::setattr(res2, "class", c("sup", class(res2)))
    }
### >2 groups
    if (ncg1 > 2) {
        ## degrees of freedom
        df1 <- seq.int(ncg1 - 1L)
        ## Subset - reduce to df1 degrees of freedom
        eMP1 <- eMP1[, .SD, .SDcols=grep("eMP_", names(eMP1))]
        ## hold results
        res3 <- data.table::data.table(array(0, dim=c(ncol(wt1), 4L)))
        data.table::setnames(res3, c("W", "chiSq", "df", "pChisq"))
        data.table::set(res3, j=1L, value=n1)
        ## We save results below as these are also used by
        ##  the tests for trend
        ## Take each column of weight, multiply by each column of eMP
        ##  then get sum of these values
        eMP1w <- apply(wt1, MARGIN=2,
                       FUN=function(wt)
                           colSums(sweep(eMP1, MARGIN=1, STATS=wt, FUN="*")))
        ## Same as above but use weight^2; then 
        cov1w <- apply(wt1, MARGIN=2,
                       FUN=function(wt)
                           rowSums(sweep(cov1, MARGIN=3, STATS=wt^2, FUN="*"),
                                   dims=2))
        dim(cov1w) <- c(ncg1, ncg1, ncol(cov1w))
        ## Subset - reduce to df1 degrees of freedom before solving
        cov1ws <- cov1w[df1, df1, ]
        cov1ws <- apply(cov1ws, MARGIN=3, FUN=solve)
        dim(cov1ws) <- c(max(df1), max(df1), length(n1))
        eMP1ss <- eMP1w[df1, ]
        ## only need subset for this calculation
        data.table::set(res3, j="chiSq",
                        value=
                            sapply(
                                seq.int(length(n1)),
                                function(i)
                                    eMP1ss[, i] %*% cov1ws[, , i] %*% eMP1ss[, i]))
        ## results
        res3[, "df" := max(df1)]
        res3[, "pChisq" := 1 - stats::pchisq(chiSq, df)]
        data.table::setattr(res3, "class", c("lrt", class(res3)))
### Tests for trend 
        ## scores - all combinations
        sAC1 <- as.matrix(expand.grid(scores, scores))
        ## scores - product of all combinations
        scoProd1 <- apply(sAC1, MARGIN=1, FUN=prod)
        data.table::set(res1, j="Q",
                        value=colSums(eMP1w * scores))
        data.table::set(res1, j="Var",
                        value=abs(apply(cov1w * scoProd1, MARGIN=3, sum)))
    }
    ## results
    res1[, "Z" := Q / sqrt(Var)]
    res1[, "pNorm" := 2 * (1 - stats::pnorm(abs(Z)))]
    res1[, "chiSq" := Q^2 / Var]
    res1[, "df" := 1]
    res1[, "pChisq" := 1 - stats::pchisq(chiSq, df)]
    data.table::setattr(res1, "class", c("lrt", class(res1)))
    ## add column for times
    ## these are unique times with at least one event
    data.table::set(wt1, j="t", value=t1)
    data.table::setattr(x, "lrw", wt1)
    if (ncg1==2) {
        data.table::setattr(x, "lrt", res1)
        data.table::setattr(x, "sup", res2)
    } else {
        data.table::setattr(x, "lrt", res3)
        res1 <- list("tft"=res1, scores=scores)
        data.table::setattr(x, "tft", res1)
    }
    print(attr(x, "lrt"))
    print(if (!is.null(attr(x, "sup"))) attr(x, "sup") else attr(x, "tft"))
    return(invisible())
}
### Helper functions
## Probability of supremum of absolute value of
##  standard Brownian motion process B(t)
## For k, 1e4 is good enough for all practical purposes
probSupBr <- function(x) {
    k <- c(0L, seq.int(1e4))
    1 - (4 / pi) * (sum(((( - 1)^k) / (2 * k + 1)) * exp(-(((pi^2) * (2 * k + 1)^2) / (8 * x^2)))))
}
## for R CMD check
S1 <- Var <- maxAbsZ <- chiSq <- Z <- df <- NULL
