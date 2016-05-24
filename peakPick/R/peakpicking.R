# Smooth peak picking -----------------------------------------------------

#' Eliminates peak pairs too close to each other.
#'
#' @param vec vector of values
#' @param candidates boolean vector representing candidate peaks
#' @param neighlim integer limit of how far apart peaks must be
#' @return boolean vector representing candidate peaks with smaller peaks of all
#'   pairs within neighlim removed
#' @keywords internal
#'
keepmax <- function (vec, candidates, neighlim) {
  repeat {
    pos <- which(candidates)
    p2pdists <- diff(pos)
    ## keep going until there is no close pair left
    if (all(p2pdists > neighlim))
      return(candidates)
    ## index smallest pairwise distance
    smallestdist <- which.min(p2pdists)
    ## pair corresponds to these positions
    pospair <- pos[smallestdist:(smallestdist + 1)]
    ## extract vector minimum for pair and map back to position
    toremove <- pospair[which.min(vec[pospair])]
    ## set that position to false (eliminate it)
    candidates[toremove] <- FALSE
  }
}


#' helper function for small peak elimination
#'
#' @param thepos integer position of peak in vector vec
#' @param vec vector of values for peakdetection
#' @param nsd numeric minimum number of standard deviations for a peak to rise
#'   above the mean of its immediate vicinity in order to be considered for a
#'   peak call
#' @param npos integer value, peak size will be estimated plus/minus npos
#'   positions from peak
#' @return boolean TRUE if the peak at position thepos is to be deleted
#' @keywords internal
#' @importFrom stats sd
#'
helperpeak <- function (thepos, vec, nsd, npos) {
  postocheck <- (thepos-npos):(thepos+npos)
  limit <- mean(vec[postocheck], na.rm=TRUE) +
    nsd*sd(vec[postocheck], na.rm=TRUE) /
    sqrt(length(which(!is.na(vec[postocheck]))))
  ifelse(vec[thepos] > limit, FALSE, TRUE)
}


#' delete small peaks from vector with peak candidates
#'
#' @param vec vector of values for peak detection
#' @param candidates boolean vector of peak candidates
#' @param nsd numeric minimum number of standard deviations for a peak to rise
#'   above the mean of its immediate vicinity in order to be considered for a
#'   peak call
#' @param npos integer value, peak size will be estimated plus/minus npos
#'   positions from peak
#' @return boolean vector of peak candidates with small peaks eliminated
#' @keywords internal
#' @importFrom utils head tail
#'
smallpeaks <- function (vec, candidates, nsd, npos) {
  ## eliminate candidates for which mean and SEM estimations are not possible
  candidates[c(head(seq_along(vec), npos), tail(seq_along(vec), npos))] <- FALSE
  pos <- which(candidates)
  if (length(pos) > 0) {
    ## go over candidate positions, return boolean for each one to delete
    todelete <- which(sapply(pos, helperpeak, vec, nsd, npos))
    ## delete detected small peaks from candidates
    candidates[pos[todelete]] <- FALSE
  }
  return(candidates)
}

#' Identifies smooth peaks in series of numbers.
#'
#' This algorithm detects peaks, smooth bumps in series of numbers.  This
#' algorithm should not be used for series containing brief spikes.  Consider
#' filtering/smoothing your data before using this algorithm.  Please refer to
#' the paper by Weber et al. for more details.
#'
#' @param mat matrix of series with series organized columnwise
#' @param neighlim integer limit for how far apart peaks must be.  Peak pairs
#'   closer than or equal to neighlim to each other have the lesser peak
#'   eliminated.
#' @param deriv.lim numeric upper limit for the estimatied derivative for a
#'   point to be considered for a peak call
#' @param peak.min.sd numeric minimum number of standard deviations for a peak
#'   to rise above the mean of its immediate vicinity in order to be considered
#'   for a peak call
#' @param peak.npos integer peak standard deviations and means will be estimated
#'   plus/minus npos positions from peak
#' @param mc.cores the number of cores to perform this computation
#' @return boolean matrix with dimensions of mat representing peaks
#' @references Weber, C.M., Ramachandran, S., and Henikoff, S. (2014).
#'   Nucleosomes are context-specific, H2A.Z-modulated barriers to RNA
#'   polymerase. Molecular Cell 53, 819-830.
#' @export
#'
peakpick <- function (mat, neighlim, deriv.lim=0.04, peak.min.sd=0.5,
                      peak.npos=10L, mc.cores=1) {
  ## coerce mat to matrix if vector
  if (is.vector(mat)) mat <- as.matrix(mat)
  ## derivative; center differences, augment with NAs at beginning and end to
  ## match original data
  der <- rbind(NA, diff(mat, lag=2)/2, NA)
  ## derivative changes excluding cases where both derivatives are exactly zero
  pos2neg <- der[1:(nrow(der) - 1), , drop=FALSE] >= 0 &
    der[2:(nrow(der)), , drop=FALSE] <= 0 &
    !(der[1:(nrow(der) - 1), , drop=FALSE] == 0 &
        der[2:(nrow(der)), , drop=FALSE] == 0)
  ## changes derivative 2->3, then the maximum can be either at 2 or 3
  ## 12 23 34 45 NA
  ## NA 12 23 34 45
  ## pos to neg either to or from point
  ## need to consider both
  signchange <- rbind(pos2neg, NA) | rbind(NA, pos2neg)
  ## additional filter; derivative magnitude
  candidates <- abs(der) < deriv.lim & signchange
  candidates[is.na(candidates)] <- FALSE
  ## eliminate small peaks
  candidates.nosmall <- parallel::mcmapply(smallpeaks, split(mat, col(mat)),
                                           split(candidates, col(candidates)),
                                           MoreArgs=list(nsd=peak.min.sd,
                                                         npos=peak.npos),
                                           SIMPLIFY=TRUE, USE.NAMES=FALSE,
                                           mc.cores=mc.cores)
  ## eliminate close peaks
  candidates.final <- parallel::mcmapply(keepmax, split(mat, col(mat)),
                               split(candidates.nosmall,
                                     col(candidates.nosmall)),
                               MoreArgs=list(neighlim=neighlim),
                               SIMPLIFY=TRUE, USE.NAMES=FALSE,
                               mc.cores=mc.cores)
  return(candidates.final)
}



# Spike detection above background ----------------------------------------

#' Helper functions for computing scaled peaks
#'
#' @param mat matrix of series with series organized columnwise
#' @param positions vector of positions to exclude from threshold computation.
#'   Internal use only; follows R's rules of matrix indexing by vector.
#' @param nsd numeric number of standard deviations for limits (see Value)
#' @return vector limits defined as means + nsd SEMs computed for the columns of
#'   mat, excluding positions from the calculation.
#' @keywords internal
#'
scaledrow <- function (mat, positions, nsd) {
  # assume matrix with vector of positions to exclude when computing sems.
  # calculate values scaled by sems for center row.
  # center determined by winlen + 1.
  matexclude <- mat
  matexclude[positions] <- NA
  means <- colMeans(matexclude, na.rm=TRUE)
  sds <- matrixStats::colSds(matexclude, na.rm=TRUE)
  # scale by sqrt(n)
  sems <- sds/sqrt(matrixStats::colCounts(!is.na(matexclude)))
  return(means + nsd*sems)
}

#' Detects spikes in series of numbers.
#'
#' This algorithm detects spikes rising above a user-specified number of
#' standard deviations numbers in a certain window.  Use this algorithm to
#' detect short spikes rather than smooth bumps in series of numbers.  Please
#' refer to the paper by Weber et al. for more details.
#'
#' @param mat matrix of series with series organized columnwise.  The algorithm
#'   treats each column separately.
#' @param roi vector of two integers (c(min, max)) defining positions in all
#'   series (rows in mat) to consider for spike detection, used together with
#'   winlen.  Must lie within the interval [2, nrow(mat) - 1].  Will be coerced
#'   to integers.
#' @param winlen integer defining the window of positions to consider for mean
#'   and sem estimation for each series.  Each estimation limits itself to the
#'   position and a plus/minus winlen positions large window.  Thus, winlen must
#'   not be chosen larger than that the windows fit within mat, given the roi.
#'   I.e. roi[1] - winlen >=1 AND roi[length(roi)] + winlen <= nrow(mat).  Will
#'   be coerced to an integer.
#' @param spike.min.sd numeric minimum number of standard deviations for a spike
#'   to rise above the mean in order to be considered for a spike call and to be
#'   excluded from the mean estimation of each subsequent iteration of the spike
#'   calling algorithm
#' @param mc.cores the number of cores do perform this calculation
#' @param verbose Boolean indicating the number of new peaks detected with each
#'   iteration.  The algorithm stops as soon as this number does not sink
#'   anymore.  Turn this on if running into problems.
#' @return boolean matrix corresponding to mat, representing spike positions.
#' @references Weber, C.M., Ramachandran, S., and Henikoff, S. (2014).
#'   Nucleosomes are context-specific, H2A.Z-modulated barriers to RNA
#'   polymerase. Molecular Cell 53, 819-830.
#' @export
detect.spikes <- function (mat, roi, winlen, spike.min.sd=3,
                           mc.cores=1, verbose=FALSE) {

  ## check matrix validity
  if (!is.matrix(mat) || nrow(mat) < 3L)
    stop("mat must be a matrix with at least 3 rows")

  ## check ROI validity
  if (!is.integer(winlen))
    winlen <- as.integer(winlen[1])
  if (!is.integer(roi))
    roi <- as.integer(roi)
  if (length(roi) != 2L || diff(roi) < 1L)
    stop("roi must contain exatly two elements, c(min, max), with max > min")
  if (roi[1] < 2L || roi[2] > nrow(mat) - 1L)
    stop("roi must lie within the interval [2, nrow(mat) - 1]")
  if (roi[1] - winlen < 1L || roi[2] + winlen > nrow(mat))
    stop("winlen too large for the roi")

  ## extract ROI
  roivec <- roi[1]:roi[2]
  roimat <- mat[roivec, ]
  ## build logical matrix, same dimensions as mat, will record spikes
  flaggedstall <- matrix(FALSE, nrow(mat), ncol(mat))
  colnames(flaggedstall) <- colnames(mat)

  ## iterative algorithm:
  changed <- integer(0)
  repeat {
    ## successively compute for all rows in roi a threshold mean + 3*sem.
    ## Exclude already detected peaks from mean and sem estimation (second
    ## argument to scaledrow)
    test <- do.call(rbind, parallel::mclapply(roivec, function (rownum) {
      scaledrow(mat[(rownum-winlen):(rownum+winlen), ],
                which(flaggedstall[(rownum-winlen):(rownum+winlen), ]),
                nsd=spike.min.sd)
    }, mc.cores=mc.cores))

    ## new iteration of spike hits
    newflaggedstall <- flaggedstall
    newflaggedstall[roivec, ] <- roimat > test

    ## if new iteration same as previous, stop
    if (identical(newflaggedstall, flaggedstall))
      return(flaggedstall)

    ## record how many new peaks detected
    changed <- append(changed, length(which(newflaggedstall != flaggedstall)))
    if (verbose)
      print(changed[length(changed)])

    ## if number of new peaks equal or more than last iteration, we reached a
    ## local minimum; stop
    if (length(changed) > 1 && diff(changed)[length(changed) - 1] >= 0)
      return(flaggedstall)

    ## otherwise, store result of iteration and start over
    flaggedstall <- newflaggedstall
  }
}
