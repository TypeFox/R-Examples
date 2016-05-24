#' Build Bootstrap Confidence Intervals for \eqn{\hat{p_k}}{p^_k}
#'
#' The function will build bootstrap confidence intervals for the bootstrap
#' estimate of  and \eqn{\mu} with a lower-bound of
#' \code{0.025} and an upper-bound of \code{0.975}.
#'
#' @inheritParams BparametersEst
#' @param bootstrap_mean A Boolean option to return the bootstrap confidence
#'  interval for the mean.
#' @param lower_bound The lower quantile for the bootstrap confidence intervals.
#' @param upper_bound The upper quantile for the bootstrap confidence intervals.
#' @references Efron, B. (1979). Bootstrap methods: another look at the
#'  jackknife. The annals of Statistics, 1-26.
#' @references Thompson, M. E., Ramirez Ramirez, L. L., Lyubchich, V. and
#'  Gel, Y. R. (2015), Using the bootstrap for statistical inference
#'  on random graphs. Can J Statistics. doi: 10.1002/cjs.11271
#' @return A list of two elements
#'  \item{p_k_CI}{This a list of length \code{length(outBootdeg$num.sam)}, one
#'    element per LSMI. Each element contains three sets of bootstrap confidence
#'    intervals for \eqn{\hat{p}_k^*} corresponding to the three estimation
#'    methods. See \code{\link{bootdeg}} for more on the three estimation methods.}
#'  \item{mean_CI}{This a list of length \code{length(outBootdeg$num.sam)}, one
#'    element per LSMI. Each element contains three sets of bootstrap confidence
#'    intervals for \eqn{\hat{\mu}} corresponding to the three estimation
#'    methods. See \code{\link{bootdeg}} for more on the three estimation methods.}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' sam.out <- Oempdegreedistrib(net = net, n.seeds = 40, n.neigh = 1, num.sam = 1)
#' outBootdeg <- bootdeg(sam.out = sam.out, n.boot = 50)
#' a <- bootCI(outBootdeg)


bootCI <- function(outBootdeg, bootstrap_mean=T, lower_bound = 0.025, upper_bound = 0.975){
  # inception apply
  # take each "list" in empd, apply CI function to each distribution in "list".
  p_k_CI <-   lapply(outBootdeg$empd,
                     FUN <- function(x) {
                       lapply(x,
                              FUN <- function(df)
                                apply(df, 2, stats::quantile, c(lower_bound, upper_bound)))})
  mean_CI <- lapply(outBootdeg$empd,
                    FUN <- function(x) {
                      lapply(x,
                             FUN <- function(df) {
                               stats::quantile(bootmeans_from_bootdegdistrib(df),
                                        c(lower_bound, upper_bound))
                               })
                      })


  if(bootstrap_mean){
    list(p_k_CI=p_k_CI, mean_CI=mean_CI)
  } else list(p_k_CI=p_k_CI)
}
