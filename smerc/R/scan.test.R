#' Spatial Scan Test
#' 
#' \code{scan.test} performs the spatial scan test of Kulldorf (1997).
#' 
#' The test is performed using the spatial scan test based on the Poisson test statistic and a fixed number of cases.  The windows are circular and extend from the observed data locations.  The clusters returned are non-overlapping, ordered from most significant to least significant.  The first cluster is the most likely to be a cluster.  If no significant clusters are found, then the most likely cluster is returned (along with a warning).
#' 
#' @param coords An \eqn{n \times 2} matrix of centroid coordinates for the regions.
#' @param cases The number of cases in each region.
#' @param pop The population size of each region.
#' @param ex The expected number of cases for each region.  The default is calculated under the constant risk hypothesis.  
#' @param type The type of scan statistic to implement.  Default is "poisson".
#' @param nsim The number of simulations from which to compute p-value.
#' @param nreport The frequency with which to report simulation progress.  The default is \code{nsim+ 1}, meaning no progress will be displayed.
#' @param ubpop The upperbound of the proportion of the total population to consider for a cluster.
#' @param alpha The significance level to determine whether a cluster is signficant.  Default is 0.05.
#' @param lonlat If lonlat is TRUE, then the great circle distance is used to calculate the intercentroid distance.  The default is FALSE, which specifies that Euclidean distance should be used.
#' @param parallel A logical indicating whether the test should be parallelized using the \code{parallel::mclapply function}.  Default is TRUE.  If TRUE, no progress will be reported.
#'
#' @return Returns a list of length two of class scan. The first element (clusters) is a list containing the significant, non-ovlappering clusters, and has the the following components:
#' \item{locids}{The location ids of regions in a significant cluster.} 
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
#' @references Waller, L.A. and Gotway, C.A. (2005).  Applied Spatial Statistics for Public Health Data.  Hoboken, NJ: Wiley.  Kulldorff, M. (1997) A spatial scan statistic. Communications in Statistics -- Theory and Methods 26, 1481-1496.
#' @examples 
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = scan.test(coords = coords, cases = floor(nydf$cases), 
#'                 pop = nydf$pop, nsim = 49, 
#'                 alpha = 0.12, lonlat = TRUE)
#' ## plot output for new york state
#' # specify desired argument values
#' mapargs = list(database = "state", region = "new york", 
#' xlim = range(out$coords[,1]), ylim = range(out$coords[,2]))
#' # needed for "state" database (unless you execute library(maps))
#' data(stateMapEnv, package = "maps") 
#' plot(out, usemap = TRUE, mapargs = mapargs)
#' 
#' # a second example to match the results of Waller and Gotway (2005)
#' # in chapter 7 of their book (pp. 220-221).
#' # Note that the 'longitude' and 'latitude' used by them has 
#' # been switched.  When giving their input to SatScan, the coords
#' # were given in the order 'longitude' and 'latitude'.
#' # However, the SatScan program takes coordinates in the order 
#' # 'latitude' and 'longitude', so the results are slightly different
#' # from the example above.
#' coords = with(nydf, cbind(y, x))
#' out2 = scan.test(coords = coords, cases = floor(nydf$cases), 
#'                   pop = nydf$pop, nsim = 49, 
#'                   alpha = 0.5, lonlat = TRUE)
#' # the cases observed for the clusters in Waller and Gotway: 117, 47, 44
#' # the second set of results match
#' c(out2$clusters[[1]]$cases, out2$clusters[[2]]$cases, out2$clusters[[3]]$cases)
scan.test = function (coords, cases, pop, ex = sum(cases)/sum(pop)*pop, 
                        type = "poisson",
                        nsim = 499, alpha = 0.1, nreport = nsim + 1, 
                        ubpop = 0.5, lonlat = FALSE, parallel = TRUE) 
{
  # argument checking
  arg_check_scan_test(coords, cases, pop, ex, nsim, alpha, 
                      nreport, ubpop, lonlat, parallel, 
                      k = 1, w = diag(nrow(coords)))
  
  # convert to propert format
  coords = as.matrix(coords)
  N = nrow(coords)
  # short names
  y = cases; e = ex
  
  if(lonlat)
  {
    d = fields::rdist.earth(coords, coords, miles = FALSE)
  }else
  {
    d = SpatialTools::dist1(coords)
  }
  
  # for each region, determine sorted nearest neighbors
  # subject to population constraint
  mynn = nnpop(d, pop, ubpop)

  # mynn = spdep::knearneigh(coords, k = (k - 1), longlat = lonlat)$nn
  # mynn = cbind(1:N, mynn)
  
  
  # display sims completed, if appropriate
  if (nreport <= nsim && !parallel) cat("sims completed: ")
  
  # determine the expected cases in/out each successive 
  # window, total number of cases, total population
  ein = unlist(lapply(mynn, function(x) cumsum(e[x])))
  ty = sum(y) # sum of all cases
  eout = ty - ein # counts expected outside the window
  popin = unlist(lapply(mynn, function(x) cumsum(pop[x])))
  
  # determine which call for simulations
  fcall = lapply
  if (parallel) fcall = parallel::mclapply
  # setup list for call
  fcall_list = list(X = as.list(1:nsim), FUN = function(i){
    # simulate new data set
    ysim = stats::rmultinom(1, size = ty, prob = e)
    # cumulate the number of cases inside the successive windows
    yin = unlist(lapply(mynn, function(x) cumsum(ysim[x])))
    # calculate all test statistics
    tall = scan.stat(yin, ein, eout, ty, type)
    # update progress
    if ((i%%nreport) == 0) cat(paste(i, ""))
    # return max of statistics for simulation
    return(max(tall))
  })
  
  # use mclapply or lapply to find max statistics for each simulation
  tsim = unlist(do.call(fcall, fcall_list), use.names = FALSE)
  
  # number of nns for each observation
  nnn = unlist(lapply(mynn, length))
  
  # factors related to number of neighbors
  # each event location possesses
  fac = rep(1:N, times = nnn)
  
  # determine yin and yout for all windows for observed data
  yin = unlist(lapply(mynn, function(x) cumsum(y[x])))
  yout = ty - yin
  
  ### calculate scan statistics for observed data
  # of distance from observation centroid
  tobs = scan.stat(yin, ein, eout, ty, type = type)
  # max scan statistic over all windows
  tscan = max(tobs)
  # observed test statistics, split by centroid in order of successive windows
  tobs_split = split(tobs, fac)
  
  # position of most likely cluster centered at each centroid
  tmax_pos = lapply(tobs_split, which.max)
  
  # determine the farthest neighbor of each centroid for
  # which the maximum scan statistic occurs for that centroid
  max_nn = mapply(function(a, b) a[b], a = mynn, b = tmax_pos)
  # index of the the windows where the max statistic occurs for each centroid
  # in the vector of all scan statistics
  tmax_idx = cumsum(c(0, nnn[-N])) + unlist(tmax_pos)
  
  # value of statistic for most likely cluster centered at each centroid
  tmax = lapply(tobs_split, max, na.rm = TRUE)
  
  # p-values associated with these max statistics for each centroid
  pvalue = unname(sapply(tmax, function(x) (sum(tsim >= x) + 1)/(nsim + 1)))
  
  # determine which potential clusters are significant
  sigc = which(pvalue <= alpha, useNames = FALSE)
  
  # if there are no significant clusters, return most likely cluster
  if(length(sigc) == 0)
  {
    sigc = which.max(tmax)
    warning("No significant clusters.  Returning most likely cluster.")
  }
  
  # which statistics are significant
  sig_tscan = unlist(tmax, use.names = FALSE)[sigc]
  # order statistics from smallest to largest
  o_sig = order(sig_tscan, decreasing = TRUE)
  # idx of significant clusters in order of significance
  sigc = sigc[o_sig]
  
  # determine the location ids in each significant cluster
  sig_regions = mapply(function(a, b) mynn[[a]][1:b], a = sigc, b = tmax_pos[sigc], SIMPLIFY = FALSE) 
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
  sig_tstat = tmax[usigc]
  sig_p = pvalue[usigc]
  sig_coords = coords[usigc,, drop = FALSE]
  sig_r = d[cbind(usigc, max_nn[usigc])]
  sig_yin = (yin[tmax_idx])[usigc]
  sig_ein = (ein[tmax_idx])[usigc]
  sig_popin = (popin[tmax_idx])[usigc]
  sig_smr = sig_yin/sig_ein
  sig_rr = (sig_yin/sig_popin)/((ty - sig_yin)/(sum(pop) - sig_popin))
  
  # reformat output for return
  clusters = vector("list", length(u))
  for(i in seq_along(clusters))
  {
    clusters[[i]]$locids = sig_regions[[i]]
    clusters[[i]]$coords = sig_coords[i,, drop = FALSE]
    clusters[[i]]$r = sig_r[i]
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

# argument checking for all scan tests
arg_check_scan_test = 
  function(coords, cases, pop, ex, nsim, alpha, nreport,
           ubpop, lonlat, parallel, k, w)
{
    if(ncol(coords) != 2) stop("coords must have two columns")
    N = nrow(coords)
    if(length(cases) != N) stop("length(cases) != nrow(coords)")
    if(!is.numeric(cases)) stop("cases should be a numeric vector")
    if(length(pop) != N) stop("length(pop) != nrow(coords)")
    if(!is.numeric(pop)) stop("pop should be a numeric vector")
    if(length(ex) != N) stop("length(ex) != nrow(coords)")
    if(!is.numeric(ex)) stop("ex should be a numeric vector")
    if(length(alpha) != 1 || !is.numeric(alpha)) stop("alpha should be a numeric vector of length 1")
    if(alpha < 0 || alpha > 1) stop("alpha should be a value between 0 and 1")
    if(length(nsim) != 1 || !is.numeric(nsim)) stop("nsim should be a vector of length 1")
    if(nsim < 1) stop("nsim should be an integer of at least 1")
    if(length(ubpop) != 1 || !is.numeric(ubpop)) stop("ubpop should be a numeric vector of length 1")
    if(ubpop<= 0 || ubpop > 1) stop("ubpop should be a value between 0 and 1")
    if(length(lonlat) != 1) stop("length(lonlat) != 1")
    if(!is.logical(lonlat)) stop("lonlat should be a logical value")
    if(length(parallel) != 1) stop("length(parallel) != 1")
    if(!is.logical(parallel)) stop("parallel should be a logical value")
    if(length(k) != 1) stop("k must have length 1")
    if(k < 1) stop("k must be an integer >= 1")
    if(!is.matrix(w)) stop("w must be a matrix")
    if(nrow(w) != ncol(w)) stop("w much be a square matrix")
    if(!is.numeric(w)) stop("w must be a numeric matrix")
    if(nrow(w) != nrow(coords)) stop("nrow(w) != nrow(coords)")
    if(floor(k) > nrow(coords)) stop("k cannot be more than the number of regions.")
}


