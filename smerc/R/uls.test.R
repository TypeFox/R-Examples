#' Upper Level Set Spatial Scan Test
#' 
#' \code{uls.test} performs Upper Level Set (ULS) spatian scan test of Patil and Taillie (2004).
#' 
#' The test is performed using the spatial scan test based on the Poisson test statistic and a fixed number of cases.  The windows are based on the Upper Level Sets proposed by Patil and Taillie (2004).  The clusters returned are non-overlapping, ordered from most significant to least significant.  The first cluster is the most likely to be a cluster.  If no significant clusters are found, then the most likely cluster is returned (along with a warning).
#' 
#' @param coords An \eqn{n \times 2} matrix of centroid coordinates for the regions.
#' @param cases The number of cases in each region.
#' @param pop The population size of each region.
#' @param w The binary spatial adjacency matrix.
#' @param ex The expected number of cases for each region.  The default is calculated under the constant risk hypothesis.  
#' @param nsim The number of simulations from which to compute p-value.
#' @param nreport The frequency with which to report simulation progress.  The default is \code{nsim+ 1}, meaning no progress will be displayed.
#' @param ubpop The upperbound of the proportion of the total population to consider for a cluster.
#' @param alpha The significance level to determine whether a cluster is signficant.  Default is 0.05.
#' @param lonlat If lonlat is TRUE, then the great circle distance is used to calculate the intercentroid distance.  The default is FALSE, which specifies that Euclidean distance should be used.
#' @param parallel A logical indicating whether the test should be parallelized using the \code{parallel::mclapply function}.  Default is TRUE.  If TRUE, no progress will be reported.
#'
#' @return Returns a list of length two of class scan. The first element (clusters) is a list containing the significant, non-ovlappering clusters, and has the the following components: 
#' \item{locids}{The location ids of regions in a significant cluster.}
#' \item{pop}{The total population in the cluser window.}
#' \item{cases}{The observed number of cases in the cluster window.}
#' \item{expected}{The expected number of cases in the cluster window.}
#' \item{smr}{Standarized mortaility ratio (observed/expected) in the cluster window.}
#' \item{rr}{Relative risk in the cluster window.}
#' \item{loglikrat}{The loglikelihood ratio for the cluster window (i.e., the log of the test statistic).}
#' \item{pvalue}{The pvalue of the test statistic associated with the cluster window.}
#' The second element of the list is the centroid coordinates.  This is needed for plotting purposes.
#' @author Joshua French
#' @importFrom parallel mclapply
#' @importFrom fields rdist.earth
#' @importFrom smacpod noc
#' @importFrom stats rmultinom
#' @export
#' @references Waller, L.A. and Gotway, C.A. (2005).  Applied Spatial Statistics for Public Health Data.  Hoboken, NJ: Wiley.  Kulldorff, M. (1997) A spatial scan statistic. Communications in Statistics -- Theory and Methods 26, 1481-1496.
#' @examples 
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = uls.test(coords = coords, cases = floor(nydf$cases), 
#'                   pop = nydf$pop, w = nyw, 
#'                   alpha = 0.12, lonlat = TRUE,
#'                   nsim = 10, ubpop = 0.1)
#' ## plot output for new york state
#' # specify desired argument values
#' mapargs = list(database = "state", region = "new york", 
#' xlim = range(out$coords[,1]), ylim = range(out$coords[,2]))
#' # needed for "state" database (unless you execute library(maps))
#' data(stateMapEnv, package = "maps") 
#' plot(out, usemap = TRUE, mapargs = mapargs)
#' 
#' data(nypoly)
#' library(sp)
#' plot(nypoly, col = color.clusters(out))
uls.test = function (coords, cases, pop, w, ex = sum(cases)/sum(pop)*pop,  
                        nsim = 499, alpha = 0.1, nreport = nsim + 1, 
                        ubpop = 0.5, lonlat = FALSE, parallel = TRUE) 
{
  # sanity checking
  arg_check_scan_test(coords, cases, pop, ex, nsim, alpha, 
                      nreport, ubpop, lonlat, parallel, 
                      k = 1, w = w)
  
  coords = as.matrix(coords)
  N = nrow(coords)
 
  y = cases
  e = ex
  ty = sum(y) # sum of all cases
  
  # display sims completed, if appropriate
  if (nreport <= nsim && !parallel) cat("sims completed: ")
  
  fcall = lapply
  if (parallel) fcall = parallel::mclapply
  fcall_list = list(X = as.list(1:nsim), FUN = function(i){
    # simulate new data set
    ysim = stats::rmultinom(1, size = ty, prob = e)
    uz = uls.zones(ysim, pop, w, ubpop = ubpop)
    # cumulate the number of cases inside the successive windows
    yin = unlist(lapply(uz, function(x) sum(ysim[x])), use.names = FALSE)
    ein = unlist(lapply(uz, function(x) sum(e[x])), use.names = FALSE)
    # calculate all test statistics
    tall = scan.stat(yin, ein, ty - ein, ty)
    # update progress
    if((i%%nreport) == 0) cat(i, " ")
    # return max of statistics for simulation
    return(max(tall))
  })
  
  # use mclapply or lapply to find max statistics for each simulation
  tsim = unlist(do.call(fcall, fcall_list), use.names = FALSE)
  
  # determine yin and yout for all windows for observed data
  uz = uls.zones(y, pop, w, ubpop)
  yin = unlist(lapply(uz, function(x) sum(y[x])), use.names = FALSE)
  yout = ty - yin
  ein = unlist(lapply(uz, function(x) sum(e[x])), use.names = FALSE)
  eout = ty - ein
  popin = unlist(lapply(uz, function(x) sum(pop[x])), use.names = FALSE)
  
  ### calculate scan statistics for observed data
  # of distance from observation centroid
  tobs = scan.stat(yin, ein, ty - ein, ty)
  # max scan statistic over all windows
  tscan = max(tobs)
 
  # p-values associated with these max statistics for each centroid
  pvalue = unname(sapply(tobs, function(x) (sum(tsim >= x) + 1)/(nsim + 1)))
  
  # determine which potential clusters are significant
  sigc = which(pvalue <= alpha, useNames = FALSE)
  
  # if there are no significant clusters, return most likely cluster
  if(length(sigc) == 0)
  {
    sigc = which.max(tobs)
    warning("No significant clusters.  Returning most likely cluster.")
  }
  
  # which statistics are significant
  sig_tscan = tobs[sigc]
  # order statistics from smallest to largest
  o_sig = order(sig_tscan, decreasing = TRUE)
  # idx of significant clusters in order of significance
  sigc = sigc[o_sig]
  
  # determine the location ids in each significant cluster
  # sig_regions = mapply(function(a, b) mynn[[a]][1:b], a = sigc, b = tmax_pos[sigc], SIMPLIFY = FALSE) 
  sig_regions = uz[sigc]
  # determine idx of unique non-overlapping clusters
  u = smacpod::noc(sig_regions)
  # return non-overlapping clusters (in order of significance)
  sig_regions = sig_regions[u]
  # unique significant clusters (in order of significance)
  usigc = sigc[u]
  
  # for the unique, non-overlapping clusters in order of significance,
  # find the associated test statistic, p-value, centroid,
  # window radius, cases in window, expected cases in window, 
  # population in window, standarized mortality ration, relative risk,
  sig_tstat = tobs[usigc]
  sig_p = pvalue[usigc]
  sig_yin = yin[usigc]
  sig_ein = ein[usigc]
  sig_popin = popin[usigc]
  sig_smr = sig_yin/sig_ein
  sig_rr = (sig_yin/sig_popin)/((ty - sig_yin)/(sum(pop) - sig_popin))
  
  # reformat output for return
  clusters = vector("list", length(u))
  for(i in seq_along(clusters))
  {
    clusters[[i]]$locids = sig_regions[[i]]
    clusters[[i]]$pop = sig_popin[i]
    clusters[[i]]$cases = sig_yin[i]
    clusters[[i]]$expected = sig_ein[i]
    clusters[[i]]$smr = sig_smr[i]
    clusters[[i]]$rr = sig_rr[i]
    clusters[[i]]$loglikrat = sig_tstat[i]
    clusters[[i]]$pvalue = sig_p[i]
  }
  outlist = list(clusters = clusters, coords = coords)
  class(outlist) = "scan"
  return(outlist)
}


