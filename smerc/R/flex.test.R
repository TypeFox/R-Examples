#' Flexibly Shaped Spatial Scan Test
#' 
#' \code{flex.test} performs the flexibly shaped spatial scan test of Tango and Takahashi (2005).
#' 
#' The test is performed using the spatial scan test based on the Poisson test statistic and a fixed number of cases.  The first cluster is the most likely to be a cluster.  If no significant clusters are found, then the most likely cluster is returned (along with a warning).
#' 
#' @param coords An \eqn{n \times 2} matrix of centroid coordinates for the regions.
#' @param cases The number of cases in each region.
#' @param pop The population size of each region.
#' @param w An \eqn{n\times n} adjacenty matrix for the regions in the study area.
#' @param k An integer indicating the maximum number of regions to inclue in a potential cluster.  Default is 10
#' @param ex The expected number of cases for each region.  The default is calculated under the constant risk hypothesis.  
#' @param type The type of scan statistic to implement.  Default is "poisson".  Alternative is "bernoulli".
#' @param nsim The number of simulations from which to compute p-value.
#' @param nreport The frequency with which to report simulation progress.  The default is \code{nsim+ 1}, meaning no progress will be displayed.
#' @param alpha The significance level to determine whether a cluster is signficant.  Default is 0.05.
#' @param lonlat If lonlat is TRUE, then the great circle distance is used to calculate the intercentroid distance.  The default is FALSE, which specifies that Euclidean distance should be used.
#' @param parallel A logical indicating whether the test should be parallelized using the \code{parallel::mclapply function}.  Default is TRUE.  If TRUE, no progress will be reported.
#'
#' @return Returns a list of length two of class scan. The first element (clusters) is a list containing the significant, non-ovlappering clusters, and has the the following components: 
#' \item{coords}{The centroid of the significant clusters.}
#' \item{r}{The radius of the window of the clusters.}
#' \item{pop}{The total population in the cluser window.}
#' \item{cases}{The observed number of cases in the cluster window.}
#' \item{expected}{The expected number of cases in the cluster window.}
#' \item{smr}{Standarized mortaility ratio (observed/expected) in the cluster window.}
#' \item{rr}{Relative risk in the cluster window.}
#' \item{loglikrat}{The loglikelihood ratio for the cluster window (i.e., the log of the test statistic).}
#' \item{pvalue}{The pvalue of the test statistic associated with the cluster window.}
#' The second element of the list is the centroid coordinates.  This is needed for plotting purposes.
#' @author Joshua French
#' @importFrom SpatialTools dist1 dist2
#' @importFrom parallel mclapply
#' @importFrom fields rdist.earth
#' @importFrom smacpod noc
#' @importFrom stats rmultinom
#' @export
#' @references Tango, T., & Takahashi, K. (2005). A flexibly shaped spatial scan statistic for detecting clusters. International journal of health geographics, 4(1), 11.  Kulldorff, M. (1997) A spatial scan statistic. Communications in Statistics -- Theory and Methods 26, 1481-1496.
#' @examples 
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = flex.test(coords = coords, cases = floor(nydf$cases),
#'                 w = nyw, k = 3,  
#'                 pop = nydf$pop, nsim = 49, 
#'                 alpha = 0.12, lonlat = TRUE)
#'                 
#' data(nypoly)
#' library(sp)
#' plot(nypoly, col = color.clusters(out))

flex.test = function(coords, cases, pop, w, k = 10, ex = sum(cases)/sum(pop)*pop, 
                        type = "poisson",
                        nsim = 499, alpha = 0.1, nreport = nsim + 1, 
                        lonlat = FALSE, parallel = TRUE) 
{
  arg_check_scan_test(coords, cases, pop, ex, nsim, alpha, 
                      nreport, 0.5, lonlat, parallel, k = k, w = w)
  coords = as.matrix(coords)
  N = nrow(coords)
  y = cases
  e = ex
  zones = flex.zones(coords, w, k, lonlat)
  if (nreport <= nsim && !parallel) 
    cat("sims completed: ")
  ein = unlist(lapply(zones, function(x) sum(e[x])), use.names = FALSE)
  ty = sum(y)
  eout = ty - ein
  fcall = lapply
  if (parallel) 
    fcall = parallel::mclapply
  fcall_list = list(X = as.list(1:nsim), FUN = function(i) {
    ysim = stats::rmultinom(1, size = ty, prob = e)
    yin = unlist(lapply(zones, function(x) sum(ysim[x])), 
                 use.names = FALSE)
    tall = scan.stat(yin, ein, ty - ein, ty, type)
    if ((i%%nreport) == 0) cat(paste(i, ""))
    return(max(tall))
  })
  tsim = unlist(do.call(fcall, fcall_list), use.names = FALSE)
  yin = unlist(lapply(zones, function(x) sum(y[x])))
  tobs = scan.stat(yin, ein, ty - ein, ty, type = type)
  tscan = max(tobs)
  fac = sapply(zones, function(x) x[1])
  tmax_pos = tapply(tobs, fac, which.max)
  fac_idx = lapply(1:N, function(i) which(fac == i))
  tmax_idx = mapply(function(idx, pos) idx[pos], idx = fac_idx, 
                    pos = tmax_pos)
  tmax = tobs[tmax_idx]
  pvalue = sapply(tmax, function(x) (sum(tsim >= x) + 1)/(nsim + 
                                                            1))
  sig = which(pvalue <= alpha)
  if (length(sig) == 0) {
    sig = which.max(tmax)
    warning("No significant clusters.  Returning most likely cluster.")
  }
  sig = sig[order(tmax[sig], decreasing = TRUE)]
  u = smacpod::noc(zones[tmax_idx[sig]])
  usig = sig[u]
  usigidx = tmax_idx[usig]
  sig_regions = zones[usigidx]
  sig_tstat = tobs[usigidx]
  sig_p = pvalue[usig]
  sig_yin = yin[usigidx]
  sig_ein = ein[usigidx]
  sig_popin = sapply(sig_regions, function(x) sum(pop[x]))
  sig_smr = sig_yin/sig_ein
  sig_rr = (sig_yin/sig_popin)/((ty - sig_yin)/(sum(pop) - 
                                                  sig_popin))
  clusters = vector("list", length(u))
  for (i in seq_along(clusters)) {
    clusters[[i]]$locids = sig_regions[[i]]
    clusters[[i]]$pop = sig_popin[i]
    clusters[[i]]$cases = sig_yin[i]
    clusters[[i]]$expected = sig_ein[i]
    clusters[[i]]$smr = sig_smr[i]
    clusters[[i]]$rr = sig_rr[i]
    clusters[[i]]$loglikrat = sig_tstat[[i]]
    clusters[[i]]$pvalue = sig_p[i]
  }
  outlist = list(clusters = clusters, coords = coords)
  class(outlist) = "scan"
  return(outlist)
}

