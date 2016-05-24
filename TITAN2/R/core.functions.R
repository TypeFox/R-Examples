#' IndVal scores based on relative abundance across, and occurrence
#' frequency within, groups
#'
#' Calculates indicator value scores using the original method of
#' Dufrene and Legendre (1997) for all taxa in a site-by-taxon
#' matrix split into two groups along an environmental gradient. A
#' modified version (indvals) calculates indicator value scores
#' using a modification of Dufrene and Legendre (1997), whereby
#' relative abundance is computed as total relative abundance across
#' groups rather than as mean relative abundance.
#'
#' The modified version of the original IndVal method was developed
#' to address potential for bias in change point identification for
#' highly skewed samples. This and the function 'indvalps' are run
#' when the argument 'ivTot' in the original TITAN function call is
#' set to TRUE.  The function 'indval' is expected to be used more
#' frequently, and 'indvals' is provided largely for comparative
#' purposes.
#'
#' @param grp A binary vector (0,1) indicating group membership for
#'   partitioning of rows in 'spp' according to a ranks of an
#'   environmental variable.
#' @param spp A site-by-taxon matrix of taxa counts at each sampling
#'   location.
#' @param allscores A logical indicating whether the list of scores
#'   returned by the function should be truncated.  The default is
#'   to return just the largest IndVal (of two, one for each side of
#'   the partition) and on which side of the partition the score
#'   occurs. If allscores is TRUE, IndVals on both sides of the
#'   partition, as well as their relative abundances are also
#'   returned.
#' @return A 2-column matrix with the same nrow as the matrix 'spp'.
#'   The first column consists of a 1 (left) or 2 (right) for each
#'   taxon indicating on which group produced the largest IndVal
#'   score, and the second column contains the actual IndVal score.
#' @references Dufrene, M. and P. Legendre. 1997. Species
#'   assemblages and indicator species: the need for a flexible
#'   asymmetrical approach. Ecological Monographs 67:345-366.
#' @seealso \code{\link{indvalp}}, \code{\link{indvalps}},
#'   \code{\link{getivz}}, \code{\link{b.getivz}}
#' @author M. Baker and R. King
#' @keywords Indicator~Species IndVal TITAN
#' @name indval












#' @rdname indval
indval <- function(grp, spp, allscores = FALSE) {
  numTxa <- ncol(spp)

  i1 = which(grp == 1)
  i2 = which(grp == 2)

  if (length(i1) > 1) {
    mnabnd1 = colMeans(spp[i1, ])
    frq1 = colMeans(spp[i1, ] > 0)
  } else {
    mnabnd1 = spp[i1, ]
    frq1 = (spp[i1, ] > 0) * 1
  }

  if (length(i2) > 1) {
    mnabnd2 = colMeans(spp[i2, ])
    frq2 = colMeans(spp[i2, ] > 0)
  } else {
    mnabnd2 = spp[i2, ]
    frq2 = (spp[i2, ] > 0) * 1
  }

  m1 <- which(mnabnd1 == 0)
  m2 <- which(mnabnd2 == 0)

  relabnd1 <- mnabnd1/(mnabnd1 + mnabnd2)
  relabnd2 <- mnabnd2/(mnabnd1 + mnabnd2)

  relabnd1[m1] <- 0
  relabnd2[m2] <- 0

  iv1 <- (relabnd1 * frq1) * 100
  iv2 <- (relabnd2 * frq2) * 100

  ivmax <- pmax(iv1, iv2)
  maxcls <- rep(0, numTxa)
  max1 <- iv1 > iv2
  max2 <- iv2 > iv1
  maxcls = maxcls + (max1 * 1)
  maxcls = maxcls + (max2 * 2)

  ivout <- cbind(maxcls, ivmax)
  if (identical(allscores, TRUE)) {
    ivout <- cbind(maxcls, ivmax, iv1, iv2, relabnd1, relabnd2)
  }
  ivout
}
































#' @rdname indval
indvals <- function(grp, spp, allscores = FALSE) {
  numTxa <- ncol(spp)

  i1 = which(grp == 1)
  i2 = which(grp == 2)

  if (length(i1) > 1) {
    smabnd1 = colSums(spp[i1, ])
    frq1 = colMeans(spp[i1, ] > 0)
  } else {
    smabnd1 = spp[i1, ]
    frq1 = (spp[i1, ] > 0) * 1
  }

  if (length(i2) > 1) {
    smabnd2 = colSums(spp[i2, ])
    frq2 = colMeans(spp[i2, ] > 0)
  } else {
    smabnd2 = spp[i2, ]
    frq1 = (spp[i2, ] > 0) * 1
  }

  m1 <- which(smabnd1 == 0)
  m2 <- which(smabnd2 == 0)

  relabnd1 <- smabnd1/(smabnd1 + smabnd2)
  relabnd2 <- smabnd2/(smabnd1 + smabnd2)

  relabnd1[m1] <- 0
  relabnd2[m2] <- 0

  iv1 <- (relabnd1 * frq1) * 100
  iv2 <- (relabnd2 * frq2) * 100

  ivmax <- pmax(iv1, iv2)
  maxcls <- rep(0, numTxa)
  max1 <- iv1 > iv2
  max2 <- iv2 > iv1
  maxcls = maxcls + (max1 * 1)
  maxcls = maxcls + (max2 * 2)

  ivout <- cbind(maxcls, ivmax)
  if (identical(allscores, TRUE)) {
    ivsout <- cbind(maxcls, ivmax, iv1, iv2, relabnd1, relabnd2)
  }
  ivout
}























#' Calculate permuted IndVal scores using a group matrix
#'
#' This function performs the same calculations as 'indval' or
#' 'indvals', but does so using matrix operations on a set of binary
#' group assignments in matrix form. Function 'indvalps' calculates
#' indicator value scores using a modification of Dufrene and
#' Legendre (1997), whereby relative abundance is computed as total
#' relative abundance across groups rather than as mean relative
#' abundance.
#'
#' Although the 'indval' function is reasonably efficient for a
#' single calculation, when repeated during permutations (default =
#' 250) and again during for each bootstrap replicate (default =
#' 500), small differences in processing time can quickly become
#' unwieldy. Even with the matrix operation, the permutation
#' accounts for > 3/4 of processing time for most data sets (due to
#' the fact that it is repeated for each bootstrap replicate).
#'
#' The output matrix does not include information about the group
#' membership of indval maxima because this information is not used
#' in the permutation procedure (i.e., only the distribution of
#' IndVal magnitudes is relevant).
#'
#' Modification of the original IndVal method was developed to
#' address potential for bias in change point identification for
#' highly skewed samples. This and the function 'indvals' are run
#' when the argument 'ivTot' in the original TITAN function call is
#' TRUE.  It is expected that 'indval' and 'indvalp' will be used
#' more commonly.
#'
#' @param grpMatrix A site-by-permutation matrix of binary (0,1)
#'   assignments indicating group membership for successive
#'   randomizations of rows of 'spp' according to a ranks of an
#'   environmental variable.
#' @param spp A site-by-taxon matrix of taxa counts at each sampling
#'   location.
#' @return A matrix (ivmax) of IndVal maxima across each partition
#'   with nrow equal to the number of permutations and ncol equal to
#'   the number of taxa in 'spp'.
#' @references Dufrene, M. and P. Legendre. 1997. Species
#'   assemblages and indicator species: the need for a flexible
#'   asymmetrical approach. Ecol. Mon. 67:345-366.
#' @references Baker, ME and RS King.  2010. A new method for
#'   detecting and interpreting biodiversity and ecological
#'   community thresholds. Methods in Ecology and Evolution 1(1):
#'   25:37.
#' @seealso \code{\link{indval}}, \code{\link{indvals}},
#'   \code{\link{getivz}}, \code{\link{b.getivz}}
#' @author M. Baker and R. King
#' @keywords Indicator~Species IndVal TITAN
#' @name indvalp














#' @rdname indvalp
indvalp <- function(grpMatrix, spp) {
  grp1Matrix <- (grpMatrix == 1) * 1
  grp2Matrix <- (grpMatrix == 2) * 1

  grp1MeansOp <- (1/colSums(grp1Matrix)) * t(grp1Matrix)
  grp2MeansOp <- (1/colSums(grp2Matrix)) * t(grp2Matrix)

  mnabnd1 <- grp1MeansOp %*% spp
  mnabnd2 <- grp2MeansOp %*% spp

  sppCnts = (spp > 0)
  frq1 <- grp1MeansOp %*% sppCnts
  frq2 <- grp2MeansOp %*% sppCnts

  zeroIdx1 <- which(mnabnd1 == 0)
  zeroIdx2 <- which(mnabnd2 == 0)

  mnabndSumInv <- 1/(mnabnd1 + mnabnd2)
  relabnd1 <- mnabnd1 * mnabndSumInv
  relabnd2 <- mnabnd2 * mnabndSumInv

  relabnd1[zeroIdx1] <- 0
  relabnd2[zeroIdx2] <- 0

  iv1 <- (relabnd1 * frq1) * 100
  iv2 <- (relabnd2 * frq2) * 100

  ivmax <- pmax(iv1, iv2)
}













#' @rdname indvalp
indvalps <- function(grpMatrix, spp) {

  grp1Matrix <- (grpMatrix == 1) * 1
  grp2Matrix <- (grpMatrix == 2) * 1

  grp1MeansOp <- (1/colSums(grp1Matrix)) * t(grp1Matrix)
  grp2MeansOp <- (1/colSums(grp2Matrix)) * t(grp2Matrix)

  smabnd1 <- t(grp1Matrix) %*% spp
  smabnd2 <- t(grp2Matrix) %*% spp

  sppCnts = (spp > 0)
  frq1 <- grp1MeansOp %*% sppCnts
  frq2 <- grp2MeansOp %*% sppCnts

  zeroIdx1 <- which(smabnd1 == 0)
  zeroIdx2 <- which(smabnd2 == 0)

  abndInv <- 1/(smabnd1 + smabnd2)
  relabnd1 <- smabnd1 * abndInv
  relabnd2 <- smabnd2 * abndInv

  relabnd1[zeroIdx1] <- 0
  relabnd2[zeroIdx2] <- 0

  iv1 <- (relabnd1 * frq1) * 100
  iv2 <- (relabnd2 * frq2) * 100

  ivmax <- pmax(iv1, iv2)
}





















#' Permutation of group membership for a single candidate partition
#'
#' Randomizes group assignments for all permutations based on ranked
#' values of the environmental gradient, then calls the appropriate
#' matrix-based IndVal function.
#'
#' This function handles the randomization portion of the
#' permutation procedure, and then serves as a wrapper for 'indvalp'
#' and 'indvalps' (depending on the value of 'ivTot'), which
#' estimate change-point distributions across all permutations.  The
#' output is the same as 'indvalp' and 'indvalps' because the
#' function simply passes their products on.
#'
#' @param grp A vector of binary (0,1) assignments indicating group
#'   membership for a single partition of nrows in 'spp' according
#'   to a ranks of an environmental variable.
#' @param spp A site-by-taxon matrix of taxa counts at each sampling
#'   location.
#' @param ivTot A logical indicating whether IndVal scores should be
#'   calculated using total relative abundance or the mean relative
#'   abundace originally proposed by Dufrene and Legendre (1997).
#' @param nPerm The number of permutations to be performed.
#' @return A matrix (ivmax) of IndVal maxima with nrow equal to the
#'   number of permutations and ncol equal to the number of taxa in
#'   'spp'.
#' @references Dufrene, M. and P. Legendre. 1997. Species
#'   assemblages and indicator species: the need for a flexible
#'   asymmetrical approach. Ecol. Mon. 67:345-366.
#' @references Baker, ME and RS King.  2010. A new method for
#'   detecting and interpreting biodiversity and ecological
#'   community thresholds. Methods in Ecology and Evolution 1(1):
#'   25:37.
#' @seealso \code{\link{indvalp}}, \code{\link{indvalps}},
#'   \code{\link{getivz}}, \code{\link{b.getivz}}
#' @author M. Baker and R. King
#' @keywords permutation TITAN
permiv <- function(grp, spp, ivTot = ivTot, nPerm = 250) {

  grpRepMatrix <- matrix(grp, nrow = length(grp), ncol = nPerm)
  grpRandMatrix <- apply(grpRepMatrix, c(2), sample)

  if (ivTot) {
    indvalps(grpRandMatrix, spp)
  } else {
    indvalp(grpRandMatrix, spp)
  }
}

























#' Performs calculation of IndVal z scores from observed and
#' permuted values
#'
#' A wrapper function that calls 'indval' or 'indvals' to obtain
#' observed IndVal scores and 'permiv' to generate permuted values,
#' then calculates z scores and associated p-values.
#'
#' This function calls two subfunctions to first calculate IndVals
#' and their associated direction (i.e., a decreasing taxon [group
#' 1] is associated with the left side of any partition, whereas an
#' increasing taxon [group 2] is associated with the right) and
#' second to develop permuted values for each candidate partition.
#' Output includes the indicator direction (1 or 2), z scores,
#' IndVal scores, and associated p value (obtained by the fraction
#' of times an observed IndVal is greater than those obtained from
#' numPerm randomizations of equivalent group sizes).  For more
#' detail regarding the relative benefits and potential drawbacks of
#' using imax=T or imax=F, see Baker and King (2013).
#'
#' The alternative function 'b.getivz' is a streamlined version used
#' within the bootstrap procedure that does not print its status to
#' the console.
#'
#' @param clss A matrix of binary (0,1) group membership based on
#'   partitions of sampling sites ranked along an environmental
#'   gradient.
#' @param spp A site-by-taxon matrix of taxa counts at each sampling
#'   location.
#' @param ivTot A logical indicating whether IndVal scores should be
#'   calculated using total relative abundance or the mean relative
#'   abundace originally proposed by Dufrene and Legendre (1997).
#' @param nPerm The number of permutations to be used by 'permiv'.
#' @param numClass The number of classes used to partition samples
#'   along the enviorinmental gradient. The default is the total
#'   number of observations less two times the minimum split size
#'   ('minSplt').
#' @param imax A logical indicating whether taxon-specific change
#'   points should be determined using IndVal maxima or z-score
#'   maxima (as in TITAN v1.0).
#' @return A matrix containing four submatrices (the first from
#'   \code{[1:numTxa,]}, the second from
#'   \code{[(numTxa+1):(2*numTxa),]}, etc.), the first two of which
#'   include indicator direction and z scores. \itemize{ \item{Group
#'   Membership }{A vector for every taxon showing decreasing (1) or
#'   increasing (2) group membership at each value of 'envcls'}
#'   \item{z scores }{A vector for every taxon showing IndVal z
#'   scores at each value of 'envcls'} \item{IndVals }{A vector for
#'   every taxon showing IndVal scores at each value of 'envcls'}
#'   \item{p values }{A vector for every taxon showing IndVal p
#'   values at each value of 'envcls'} }
#' @references Baker, ME and RS King.  2010. A new method for
#'   detecting and interpreting biodiversity and ecological
#'   community thresholds. Methods in Ecology and Evolution 1(1):
#'   25:37.
#' @references Baker ME and RS King. 2013. Of TITAN and straw men:
#'   an appeal for greater understanding of community data.
#'   Freshwater Science 32(2):489-506.
#' @seealso \code{\link{indval}}, \code{\link{indvals}},
#'   \code{\link{indvalp}}, \code{\link{indvalps}},
#'   \code{\link{permiv}}, \code{\link{ivzsums}}
#' @author M. Baker and R. King
#' @keywords TITAN Threshold~Indicator~Taxa~Analysis Indicator~Value
#' @name getivz









#' @rdname getivz
getivz <- function(clss, spp, ivTot = ivTot, nPerm = 250, numClass = numClass,
  imax = imax) {
  numTxa <- ncol(spp)
  numUnit <- nrow(spp)
  print("Calculating observed IndVal maxima and class values")
  ## observed IndVal maxima and class values
  if (ivTot) {
    obs <- apply(clss, 2, indvals, spp = spp)
    print("Calculating IndVals using total relative abundance")
  } else {
    obs <- apply(clss, 2, indval, spp = spp)
    print("Calculating IndVals using mean relative abundance")
  }

  print("Permuting IndVal scores")
  ## Permuted IndVal scores
  spp = as.matrix(spp)
  piv <- array(apply(clss, 2, permiv, spp = spp, ivTot = ivTot,
    nPerm = nPerm), c(nPerm, numTxa, numClass))

  ## Calculate IndVal z-scores
  mniv <- colMeans(piv, dims = 1, na.rm = T)
  sdiv <- apply(piv, c(2, 3), sd)
  sdiv[which(sdiv < 0.1)] = 1
  obsiv <- obs[(numTxa + 1):(numTxa * 2), ]
  obsmx <- obs[1:numTxa, ]
  ivz <- (obsiv - mniv)/(sdiv)

  ## Calculate p-value
  iv.p <- matrix(NA, nrow(obsiv), ncol(obsiv))
  for (j in 1:ncol(obsiv)) {
    for (i in 1:nrow(obsiv)) {
      iv.p[i, j] <- (sum(piv[, i, j] >= obsiv[i, j], na.rm = TRUE) +
        1)/nPerm
    }
  }

  rbind(obsmx, ivz, obsiv, iv.p)
}
















#' @rdname getivz
b.getivz <- function(clss, spp, ivTot = ivTot, nPerm = nPerm, numClass = numClass,
  imax = imax) {
  numTxa <- ncol(spp)
  numUnit <- nrow(spp)
  ## observed IndVal maxima and class values
  if (ivTot) {
    obs <- apply(clss, 2, indvals, spp = spp)
  } else {
    obs <- apply(clss, 2, indval, spp = spp)
  }
  ## Permuted IndVal scores
  spp = as.matrix(spp)
  piv <- array(apply(clss, 2, permiv, spp = spp, ivTot = ivTot,
    nPerm = nPerm), c(nPerm, numTxa, numClass))

  ## Calculate IndVal z-scores
  mniv <- colMeans(piv, dims = 1, na.rm = T)
  sdiv <- apply(piv, c(2, 3), sd)
  sdiv[which(sdiv < 0.1)] = 1
  obsiv <- obs[(numTxa + 1):(numTxa * 2), ]
  obsmx <- obs[1:numTxa, ]
  ivz <- (obsiv - mniv)/(sdiv)

  ## Calculate p-value
  iv.p <- matrix(NA, nrow(obsiv), ncol(obsiv))
  for (j in 1:ncol(obsiv)) {
    for (i in 1:nrow(obsiv)) {
      iv.p[i, j] <- (sum(piv[, i, j] >= obsiv[i, j], na.rm = TRUE) +
        1)/nPerm
    }
  }

  rbind(obsmx, ivz, obsiv, iv.p)
}




















#' Sum IndVal z scores across taxa
#'
#' This function uses the output of 'getivz' (and, optionally,
#' similar output from the bootstrap procedure) to sum z scores
#' across taxa associated with each indicator direction (- or +) at
#' each level of the environmental gradient.
#'
#' The function selects taxa identified as either increasing or
#' decreasing at each level of the environmental gradient and
#' combines their IndVal z scores to generate an assemblage-wide
#' sum.  The sum(z) scores are interpreted as the magnitude of
#' community change at each level of the environmental gradient.
#'
#' The alias 'ivzsums.f' uses the output of 'getivz' to sum z scores
#' across taxa filtered by purity and reliability associated with
#' each indicator direction (- or +) at each level of the
#' environmental gradient.  Filters are provided by the final column
#' ("filter") in the 'sppmax' table that is part of each TITAN
#' object and an argument for 'ivzsums.f'.  All taxa with a value of
#' either 1 or 2 are pure and reliable decreasers or increasers,
#' respectively, and are selected for summation.  Filtered sums are
#' used by the 'plot.sumz' function to create estimates of robust
#' community change more precise than the original unfiltered sum(z)
#' in TITAN v1.0.
#'
#' @param allivz The output matrix from the function 'getivz' that
#'   contains four submatrices, the first two of which include
#'   indicator direction and z scores.
#' @param sppmax A completed summary output table for all taxa, used
#'   to ascertain which taxa should be filtered and which retained.
#' @return A matrix of two (z- and z+) parallel vectors with length
#'   (nrow) equal to the number of candidate partitions of an
#'   environmental gradient.
#' @references Baker, ME and RS King.  2010. A new method for
#'   detecting and interpreting biodiversity and ecological
#'   community thresholds. Methods in Ecology and Evolution 1(1):
#'   25:37.
#' @seealso \code{\link{getivz}}, \code{\link{plotSumz}}
#' @author M. Baker and R. King
#' @keywords TITAN sum(z)
#' @name ivzsums



#' @rdname ivzsums
ivzsums <- function(allivz) {
  numTxa = nrow(allivz)/4
  tmp1 <- allivz[1:numTxa, ] == 1
  tmp2 <- allivz[1:numTxa, ] == 2

  ivz <- allivz[(numTxa + 1):(numTxa * 2), ]

  s1 = colSums(ivz * tmp1, na.rm = T)
  s2 = colSums(ivz * tmp2, na.rm = T)

  sumivz <- cbind(s1, s2)
}



#' @rdname ivzsums
ivzsums.f <- function(allivz, sppmax) {
  numTxa = nrow(allivz)/4
  ivz <- allivz[(numTxa + 1):(numTxa * 2), ]
  tmp1 <- which(sppmax[, 16] == 1)
  tmp2 <- which(sppmax[, 16] == 2)
  s1 = NA
  s2 = NA
  if (length(tmp1) > 1) {
    s1 = colSums(ivz[tmp1, ], na.rm = T)
  } else {
    if (length(tmp1) > 0) {
      s1 = ivz[tmp1, ]
    }
  }
  if (length(tmp2) > 1) {
    s2 = colSums(ivz[tmp2, ], na.rm = T)
  } else {
    if (length(tmp2) > 0) {
      s2 = ivz[tmp2, ]
    }
  }
  sumivz <- cbind(s1, s2)
}

# print('Core Function definition complete')
