#' Screening of site by taxa matrix
#'
#' Screens site by taxa matrix for potential minimum occurrence and
#' frequency problems.
#'
#' This function provides a warning for low numbers of observations
#' and data sets with too few observations to make any reasonable
#' interpretation.  It also warns users when 100 percent occurence
#' is detected and how much of the data set meets this criterion, as
#' IndVal based analyses are not ideal for this type of data and we
#' might recommend other techniques (Baker and King 2013). The
#' function also assess whether occurrence frequencies for any taxa
#' fall below 3, and call attention to very small minimum split
#' sizes.
#'
#' @param txa A site by taxon matrix of sampled counts at each
#'   sampling location.
#' @param minSplt The minimum split size used for partitioning.  The
#'   default is to use the argument form the original TITAN function
#'   call.
#' @return A site by taxon matrix of sampled counts at each sampling
#'   location.
#' @references Baker, ME and RS King.  2010. A new method for
#'   detecting and interpreting biodiversity and ecological
#'   community thresholds.  Methods in Ecology and Evolution 1(1):
#'   25:37.
#' @references Baker ME and RS King. 2013. Of TITAN and straw men:
#'   an appeal for greater understanding of community data.
#'   Freshwater Science 32(2):489-506.
#' @author M. Baker and R. King
#' @keywords TITAN ~kwd2
txa.screen <- function(txa, minSplt = minSplt) {
  taxa <- data.matrix(txa)
  numUnit <- nrow(txa)
  if (numUnit < 10) {
    stop("Number of observations too small")
  }
  if (numUnit < 20) {
    warning("Low number of observations")
  }
  ocrnc <- colSums(taxa > 0, na.rm = TRUE)
  minTaxa <- min(ocrnc)
  if (max(ocrnc) == numUnit) {
    print(paste("100% occurrence detected ", length(which(ocrnc ==
      numUnit)), " times (", round((length(which(ocrnc ==
      numUnit))/numUnit) * 100, digits = 1), "% of taxa), use of TITAN less than ideal for this data type",
      sep = ""))
  }
  if (minTaxa < 3) {
    stop("Minimum taxon occurrence frequency is 3, your data do not meet this criterion")
  }

  if (minSplt < 5) {
    if (minSplt < 3) {
      stop("Minimum split size should be greater than 2, preferably 5 or more; please set minSplt to larger value and try again")
    } else {
      print("WARNING: minimum split size is preferably 5 or more; please double check minSplt settings")
    }
  }


  print("Taxa frequency screen complete")
  taxa
}































#' Partitions environmental gradient for analysis
#'
#' This function compares the number of records in the environmental
#' vector and the umber of rows in the taxa matrix to ensure they
#' are compatible.  It then uses the rank order of environmental
#' values at each sampling location to determine group membership of
#' all sampling sites at each level of partitioning.
#'
#' @param env A vector of values for each sampling location along
#'   the environmental gradient.
#' @param taxa A site by taxon matrix containing observed counts at
#'   each sampling location.
#' @param minSplt The  minimum split size used in binary
#'   partitioning.  The default is to use the argument from the
#'   original TITAN function call.
#' @return A list of seven objects: \itemize{ \item{env}{an
#'   environmental vector} \item{numUnit}{the number of sample units
#'   in env} \item{numTxa}{the number of distinct taxonomic units in
#'   taxa} \item{numClass}{the number of candidate partitions
#'   (numUnit-2*minSplt)} \item{srtEnv}{a sorted version of the
#'   environmental vector} \item{envcls}{a vector of environmental
#'   values used to distinguish partitions} \item{eclass}{a matrix
#'   of group membership relative to each partition in envcls} }
#' @references Baker, ME and RS King.  2010. A new method for
#'   detecting and interpreting biodiversity and ecological
#'   community thresholds.  Methods in Ecology and Evolution 1(1):
#'   25:37.
#' @author M. Baker and R. King
#' @keywords TITAN
env.part <- function(env, taxa, minSplt = minSplt) {
  nUnit <- length(env)
  numUnit <- nrow(taxa)
  if (nUnit != numUnit) {
    stop("Number of sites not equal between env vector and taxa matrix")
  }
  env <- as.matrix(env)
  rankEnv <- rank(env, ties.method = "random")
  srtEnv <- sort(env)
  srtEnv2 <- sort(rankEnv)
  numTxa = ncol(taxa)
  envcls <- srtEnv[(minSplt):(numUnit - minSplt)]
  numClass = length(envcls)
  eclass <- matrix(NA, numUnit, numClass)
  for (c in 1:numClass) {
    eclass[, c] <- ((rankEnv > srtEnv2[minSplt + (c - 1)]) *
      1) + 1
  }
  print("Determining partitions along gradient")
  e.part <- list(env, numUnit, numTxa, numClass, srtEnv, envcls,
    eclass)
}


























#' Summarizes results from TITAN's analaysis of observed data
#'
#' This function populates the first seven columns of the 'sppmax'
#' output table with results summarized for each taxon from the
#' output of 'getivz'.
#'
#' The items summarized for each taxon include (1) the env value at
#' the IndVal maximum, (2) the env value at the z-score maximum, (3)
#' the taxon's occurrence frequency, (4) group membership (decreaser
#' or increaser), (5) the observed IndVal score, (6) the probability
#' that the taxon's IndVal score could be generated by random
#' chance, (7) the z score. As the values 3-7 are computed for every
#' taxon at each candidate change point, values for the table are
#' determined by the maximum (IndVal or z score) indicated by the
#' value of 'imax'.  For further detail regarding the advantages and
#' disadvantages of using imax=T or imax=F, see Baker and King
#' (2013).
#'
#' @param ivzScores The product of the 'getivz' function.  A data
#'   matrix comprised of four submatrices including group
#'   membership, z scores, IndVals, and p values.
#' @param taxa A site by taxon matrix with counts observed at each
#'   sampling location.
#' @param srtEnv A sorted version of the environmental gradient.
#' @param minSplt The minimum split size used to determine
#'   partitions along the environmental gradient.  The defualt is to
#'   use the argument from the original TITAN function call.
#' @param imax A logical indicating whether taxon-specific change
#'   points should be determined using IndVal maxima or z-score
#'   maxima (as in TITAN v1.0).
#' @return The function output consists of a list of three objects:
#'   \itemize{ \item{obs1 }{a logical indicating decreasing taxa}
#'   \item{obs2 }{a logical indicating increasing taxa} \item{sppmax
#'   }{a partially completed summary output table for all taxa} }
#' @references Baker, ME and RS King.  2010. A new method for
#'   detecting and interpreting biodiversity and ecological
#'   community thresholds.  Methods in Ecology and Evolution 1(1):
#'   25:37.
#' @references Baker ME and RS King. 2013. Of TITAN and straw men:
#'   an appeal for greater understanding of community data.
#'   Freshwater Science 32(2):489-506.
#' @author M. Baker and R. King
#' @seealso \code{\link{getivz}}, \code{\link{titan}}
#' @keywords TITAN
obs.summ <- function(ivzScores, taxa, srtEnv, minSplt = minSplt,
  imax = imax) {

  ## Prep output table sppmax
  numTxa = ncol(taxa)
  sppmax <- matrix(NA, numTxa, 16)
  rownames(sppmax) <- colnames(taxa)
  colnames(sppmax) <- c("ienv.cp", "zenv.cp", "freq", "maxgrp",
    "IndVal", "obsiv.prob", "zscore", "5%", "10%", "50%", "90%",
    "95%", "purity", "reliability", "z.median", "filter")

  print("Summarizing Observed Results")
  ## Obtain env values at IndVal, z-score maxima
  max.ival <- apply(ivzScores[((numTxa * 2) + 1):(numTxa * 3),
    ], 1, which.max)
  max.zval <- apply(ivzScores[(numTxa + 1):(numTxa * 2), ], 1,
    which.max)

  ## Store direction of observed IndVals
  obs1 <- rep(NA, numTxa)
  obs2 <- rep(NA, numTxa)
  sppmax[, 16] = 0

  ## Write imax, zmax, freq, maxgrp, IndVal, iv.prob, zscores to
  ## table
  for (i in 1:numTxa) {
    sppmax[i, 1] <- (srtEnv[minSplt + (max.ival[i] - 1)] +
      srtEnv[minSplt + max.ival[i]])/2
    sppmax[i, 2] <- (srtEnv[minSplt + (max.zval[i] - 1)] +
      srtEnv[minSplt + max.zval[i]])/2
    sppmax[i, 3] <- sum(taxa[, i] > 0, na.rm = TRUE)
    if (imax == FALSE) {
      if (i == 1) {
        print("Estimating taxa change points using z-score maxima")
      }
      sppmax[i, 4] <- ivzScores[i, max.zval[i]]
      sppmax[i, 5] <- ivzScores[((numTxa * 2) + i), max.zval[i]]
      sppmax[i, 6] <- ivzScores[((numTxa * 3) + i), max.zval[i]]
      sppmax[i, 7] <- ivzScores[((numTxa * 1) + i), max.zval[i]]
      sppmax[, 7] <- round(sppmax[, 7], 2)
      sppmax[, 5] <- round(sppmax[, 5], 2)
      obs1[i] <- ivzScores[i, max.zval[i]] == 1
      obs2[i] <- ivzScores[i, max.zval[i]] == 2
    } else {
      if (i == 1) {
        print("Estimating taxa change points using raw IndVal maxima")
      }
      sppmax[i, 4] <- ivzScores[i, max.ival[i]]
      sppmax[i, 5] <- ivzScores[((numTxa * 2) + i), max.ival[i]]
      sppmax[i, 6] <- ivzScores[((numTxa * 3) + i), max.ival[i]]
      sppmax[i, 7] <- ivzScores[((numTxa * 1) + i), max.ival[i]]
      sppmax[, 7] <- round(sppmax[, 7], 2)
      sppmax[, 5] <- round(sppmax[, 5], 2)
      obs1[i] <- ivzScores[i, max.ival[i]] == 1
      obs2[i] <- ivzScores[i, max.ival[i]] == 2
    }
  }
  list(obs1, obs2, sppmax)
}





























#' Summarizes the results of the community-level sum(z) analysis
#'
#' This function populates the sumz.cp table using the resuls from
#' function 'ivzsums' and, if 'boot'=TRUE, calls 'ivsums.f' to
#' compute ivz sums filtered by pure and reliable taxa.
#'
#' The function locates the env values of sum(z) maxima, then if
#' 'boot'=TRUE, locates the the env value of the filtered sum(z) and
#' provides bootstrap quantiles of both filtered and unfiltered
#' distributions.
#'
#' @param ivzScores The product of the 'getivz' function.  A data
#'   matrix comprised of four submatrices including group
#'   membership, z scores, IndVals, and p values.
#' @param ivz The product of the 'ivzsums' function.  A data matrix
#'   comprised of two parallel vectors of sum(z-) and sum(z+)
#'   scores.
#' @param srtEnv A sorted version of the environmental gradient.
#' @param sppmax The taxon-specific summary output table from TITAN,
#'   passed to 'ivzsums.f'.
#' @param maxSumz A vector of sum(z) maxima across bootstrap
#'   replicates.
#' @param maxFsumz A vector of sum(z) maxima filtered by pure and
#'   reliable taxa across bootstrap replicates.
#' @param minSplt The minimum split size used to partition the
#'   environmental gradient.  The default is to use the argument
#'   specified by the original TITAN function call.
#' @param boot A logical indicating whether the bootstrap procedure
#'   should be implemented.  The default is to use the argument
#'   specified by the original TITAN function call.
#' @return A list with two objects: \itemize{ \item{sumz.cp }{A
#'   second summary output table from TITAN to accompany 'sppmax'.}
#'   \item{ivz.f }{The product of the 'ivzsums.f' function.  A data
#'   matrix comprised of two parallel vectors of filtered sum(z-)
#'   and sum(z+) scores.} }
#' @references Baker, ME and RS King.  2010. A new method for
#'   detecting and interpreting biodiversity and ecological
#'   community thresholds.  Methods in Ecology and Evolution 1(1):
#'   25:37.
#' @references King, RS and ME Baker  2010. Considerations for
#'   identifying and interpreting ecological community thresholds.
#'   Journal of the North American Benthological Association
#'   29(3):998-1008.
#' @author M. Baker and R. King
#' @seealso \code{\link{ivzsums}}, \code{\link{ivzsums.f}},
#'   \code{\link{getivz}}, \code{\link{titan}}
#' @keywords TITAN sum(z)
sumz.tab <- function(ivzScores, ivz, srtEnv, sppmax, maxSumz = maxSumz,
  maxFsumz = maxFsumz, minSplt = minSplt, boot = boot) {
  sumz.cp <- matrix(NA, 4, 6)
  colnames(sumz.cp) <- c("cp", "0.05", "0.10", "0.50", "0.90",
    "0.95")
  rownames(sumz.cp) <- c("sumz-", "sumz+", "fsumz-", "fsumz+")
  sumz.cp[1, 1] <- (srtEnv[minSplt + (which.max(ivz[, 1]) - 1)] +
    srtEnv[minSplt + which.max(ivz[, 1])])/2
  sumz.cp[2, 1] <- (srtEnv[minSplt + (which.max(ivz[, 2]) - 1)] +
    srtEnv[minSplt + which.max(ivz[, 2])])/2
  ivz.f <- 0
  if (boot) {
    ivz.f <- ivzsums.f(ivzScores, sppmax)
    if (!is.na(sum(ivz.f[, 1]))) {
      sumz.cp[3, 1] <- (srtEnv[minSplt + (which.max(ivz.f[,
        1]) - 1)] + srtEnv[minSplt + which.max(ivz.f[,
        1])])/2
    }
    if (!is.na(sum(ivz.f[, 2]))) {
      sumz.cp[4, 1] <- (srtEnv[minSplt + (which.max(ivz.f[,
        2]) - 1)] + srtEnv[minSplt + which.max(ivz.f[,
        2])])/2
    }
    sumz.cp[1, 2:6] <- quantile(maxSumz[, 1], probs = c(0.05,
      0.1, 0.5, 0.9, 0.95), na.rm = TRUE)
    sumz.cp[2, 2:6] <- quantile(maxSumz[, 2], probs = c(0.05,
      0.1, 0.5, 0.9, 0.95), na.rm = TRUE)
    sumz.cp[3, 2:6] <- quantile(maxFsumz[, 1], probs = c(0.05,
      0.1, 0.5, 0.9, 0.95), na.rm = TRUE)
    sumz.cp[4, 2:6] <- quantile(maxFsumz[, 2], probs = c(0.05,
      0.1, 0.5, 0.9, 0.95), na.rm = TRUE)
  }
  sumzTabList <- list(sumz.cp, ivz.f)
}
