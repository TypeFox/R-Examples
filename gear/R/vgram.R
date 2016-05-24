#' Sample semivariogram
#' 
#' \code{vgram} calculates the sample semivariogram.
#' 
#' Note that the directions may be different from other packages (e.g., \code{gstat} or \code{geoR} packages) because those packages calculate angles clockwise from the y-axis, which is a convention frequently seen in geostatistics (e.g., the GSLIB software library).
#' 
#' Additionally, note that calculating the sample semivariogram for the residuals of lm(response ~ 1) will produce identical results to simply computing the sample semivariogram from the original response.  However, if a trend is specified (the righthand side of ~ has non-trival covariates), then the sample semivariogram of the residuals will differ from that of the original response.  A trend should be specified when the mean is non-stationary over the spatial domain.    
#' 
#' @param formula A formula describing the relationship between the response and any covariates of interest, e.g., response ~ 1.  The variogram is computed for the residuals of the linear model \code{lm(formula, data)}.
#' @param data A data.frame, SpatialPointsDataFrame, SpatialPixelsDataFrame, or SpatialGridDataFrame object.
#' @param coords A formula specifiying (on the righthand side of ~) the coordinates in \code{data}.  Default is NULL, but this must be specified if \code{data} is of class \code{data.frame}.  This will probably look something like "~ x + y".  
#' @param nbins The number of bins (tolerance regions) to use when estimating the sample semivariogram.
#' @param maxd The maximum distance used when calculating the semivariogram.  Default is NULL, in which case half the maximum distance between coordinates is used.
#' @param angle A single value (in degrees) indicating the starting direction for a directional variogram.  The default is 0.
#' @param ndir The number of directions for which to calculate a sample semivariogram.  The default is 1, meaning calculate an omnidirection semivariogram.
#' @param type The name of the estimator to use in the estimation process.  The default is "standard", the typical method-of-moments estimator.  Other options include "cressie" for the robust Cressie-Hawkins estimator, and "cloud" for a semivariogram cloud based on the standard estimator.  If "cloud" is specified, the \code{nbins} argument is ignored.
#' @param npmin The minimum number of pairs of points to use in the semivariogram estimator.  For any bins with fewer points, the estimate for that bin is dropped.
#' @param longlat A logical indicating whether Euclidean (\code{FALSE}) or Great Circle distance (WGS84 ellipsoid) (\code{longlat = TRUE}) should be used.  Default is \code{FALSE}.
#' 
#' @return Returns a vgram object with components: 
#' @author Joshua French
#' @export
#' @importFrom parallel mclapply
#' @importFrom sp coordinates spDists
#' @importFrom stats lm
#' @examples 
#' data(co)
#' v = vgram(Al ~ 1, co, ~ easting + northing)
#' plot(v)
#' v2 = vgram(Al ~ 1, co, ~ easting + northing, angle = 22.5, ndir = 4)
#' plot(v2)

vgram = function(formula, data, coords = NULL, nbins = 10, 
                 maxd = NULL, angle = 0, ndir = 1, 
                     type = "standard", npmin = 2, 
                 longlat = FALSE)
{
  vgram_check_args(formula, data, coords, nbins, maxd, 
                   angle, ndir, type, npmin)
  id = all.vars(formula)[1]
  if(is.data.frame(data)) sp::coordinates(data) = coords
  # extract coordinates, determine number of coordinates
  coords = sp::coordinates(data)
  N = nrow(coords)
  
  # extract (residual) response
  # y = slot(data, "data")[,id]
  y = stats::lm(formula, data = data)$resid
  
  dirname = "omnidirectional"
  if(ndir >= 2) dirname = "directional"
  cat(paste("Calculating", dirname,"sample semivariogram for", id, "variable using", type, "estimator"))
  
  # create indexes to use in building of sample semivariograms
  idx1 = rep(1:(N-1), times = (N-1):1)
  idx2 = unlist(sapply(1:(N-1), function(i) (1:N)[-(1:i)]))

  # distances for unique pairs of points
  # d2 = c(dist(coords)) 
  d = sp::spDists(as.matrix(coords), longlat = longlat)
  d = c(d[lower.tri(d)])
  
  # difference of unique pairs of points
  diff = y[idx1] - y[idx2]
  
  # determine maximum distance of estimator
  if(is.null(maxd)) maxd = max(d)/2
  # determine distances of tolerance regions
  # for sample semivariogram estimator
  bindist = seq(0, maxd, len = nbins + 1)
  # determine classification cuts for unique pairs of points by distance
  dcuts = cut(d, breaks = bindist)

  # default angle_dir, bin_angles
  angle_dir = angle
  bin_angles = "[0, 180)"
  if(ndir > 1)
  {
    # determine change between angle directions
    angle_change = 180/ndir
    # determine main angles
    angle_dir = (angle + 0:(ndir-1) * angle_change)%%180 
    # determine tolerance bins for angles
    bin_angles = angle + (seq(-1, (2 * ndir - 1), by = 2))*angle_change/2
    # if any angles are negative, correct them by rotating 180 degrees
    which_neg = which(bin_angles < 0)
    bin_angles[which_neg] = (180 - bin_angles[which_neg])
    bin_angles = sort(bin_angles)
    # calculate angles between observed data
    obs_angles = angle2d(coords[idx1, ], coords[idx2, ]) %% 180
    min_angle  = min(bin_angles)
    obs_angles[obs_angles < min_angle] = obs_angles[obs_angles < min_angle] + 180 
    # determine split by angle
    acuts = cut(obs_angles, bin_angles, right = FALSE) # don't include 180
    # determine indices of split
    split_idx = split(1:length(d), acuts)
    # determine semivariograms for each direction
    semi = parallel::mclapply(split_idx, 
                              function(idx)
                              {
                                omni_semivariogram(d[idx], diff[idx], dcuts[idx], npmin, type)
                              })
    # determine angle_bin
    angle_name = rep(names(semi), times = unlist(lapply(semi, nrow), use.names = FALSE))  
    # bind semivariograms from each direction
    semi = do.call(rbind.data.frame, semi)
    # remove ugly row names
    row.names(semi) = NULL
    # add name of angle bin
    semi$angle = angle_name
  }else
  {
    semi = omni_semivariogram(d, diff, dcuts, npmin, type)
  }
  out = list(id = id, 
             nbins = nbins,
             bindist = bindist,
             maxd = maxd,
             angle = angle_dir,
             ndir = ndir,
             binangles = bin_angles,
             type = type,
             semivariogram = semi,
             object = data)
  class(out) = "vgram"
  return(out)
}

vgram_check_args = function(formula, data, coords, nbins, maxd, 
                            angle, ndir, type, npmin)
{
  if(class(formula) != "formula") stop("formula is not of class formula")
  if(is.null(formula[[2]]) || length(formula[[2]]) > 2) stop("formula should contain a single variable to the left of ~")
  if(!is.element(class(data), c("data.frame", "SpatialPointsDataFrame", "SpatialGridDataFrame", "SpatialPixelsDataFrame")))
  {
    stop("data not of appropriate class.  Should be of class data.frame, SpatialPointsDataFrame, SpatialGridDataFrame, or SpatialPixelsDataFrame.")
  }
  if(min(is.element(all.vars(formula), names(as.data.frame(data)))) == 0)
    stop("some of the variables in formula are not in data")
  if(is.data.frame(data))
  {
    if(is.null(coords)) stop("coords cannot be NULL when data is a data.frame") 
    if(class(coords) != "formula") stop("coords must be a formula")
    if(min(is.element(all.vars(coords), names(as.data.frame(data)))) == 0)
      stop("some of the variables in coords are not in data")
  }
  if(nbins < 1 || length(nbins) != 1 || !is.numeric(nbins))
    stop("nbins should be a single integer >= 1")
  if(!is.null(maxd))
  {
    if(maxd <= 0 || length(maxd) != 1 || !is.numeric(maxd))
      stop("maxd should be a single integer >= 1")
  }
  if(length(angle) != 1 || !is.numeric(angle) || angle < 0)
    stop("angle should be an single value >= 0")
  if(ndir < 1 || length(ndir) != 1 || !is.numeric(ndir))
    stop("ndir should be a single integer >= 1")
  if(length(type) != 1) stop("type should be a single name")
  if(!is.element(type, c("standard", "cressie", "cloud")))
    stop(paste(type, "is not a valid type"))
  if(npmin < 1 || length(npmin) != 1 || !is.numeric(npmin))
    stop("npmin should be a single integer >= 1")
}
    
# determine omnidireciton semivariogram based on 
# based on distances, kernel (the numerator of the estimator),
# the cuts (tolerance regions), and minimum numer of pairs allowed
omni_semivariogram = function(d, diff, cuts, npmin, type)
{
  np = unlist(lapply(split(d, cuts), length), use.names = FALSE)
  
  if(type == "standard")
  {
    semivariance = unlist(lapply(split(diff^2, cuts), mean), use.names = FALSE)/2  
  }else if(type == "cressie")
  {
    semivariance = unlist(lapply(split(sqrt(abs(diff)), cuts), mean), use.names = FALSE)^4/2/(0.457 + 0.494/np)
  }else if(type == "cloud")
  {
    return(data.frame(distance = d, semivariance = diff^2/2, np = 1))
  }
  
  distance  = unlist(lapply(split(d, cuts), mean), use.names = FALSE)
  which_small = which(np < npmin)
  if(length(which_small) > 0)
  {
    return(data.frame(distance = distance[-which_small], 
                semivariance = semivariance[-which_small],
                np = np[-which_small]))
  }else
  {
    return(data.frame(distance = distance, 
                semivariance = semivariance,
                np = np))     
  } 
}