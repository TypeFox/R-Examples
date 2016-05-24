#' Calculate bootstrapped IndVal z scores
#'
#' This function implements resampling (with replacement) of the
#' observed environmental gradient and site-by-taxon matrix, and
#' then calls the function \code{b.getivz} to obtain
#' bootstrapped scores.
#'
#' Four pieces of information are obtained from every taxon during
#' each bootstrap replicate. If the argument 'imax' is TRUE,
#' bootstrapped change points are identified based on IndVal maxima,
#' whereas if 'imax' is FALSE z-score maxima are used instead.  In
#' addition to the IndVal or z score maxima, the value of the
#' environmental gradient, the indicator direction, and the p value
#' are also retained for that point.
#'
#' In addition to the above metric matrix for each taxon (1), the z
#' scores across all candidate chnage points are retained from each
#' replicate (2), as well as the response direction (maxgrp) (3) and
#' a sorted version of resampled environmental values (4).  These
#' four items are combined as a list object.
#'
#' @param bSeq An index used to determine the sequence number of the
#'   current bootstrap replicate.
#' @param env An environmental gradient.
#' @param taxa A site-by-taxon matrix of taxa counts at each
#'   sampling location.
#' @param ivTot A logical indicating whether IndVal scores should be
#'   calculated using total relative abundance or the mean relative
#'   abundace originally proposed by Dufrene and Legendre (1997).
#'   The default is to pass on the argument from the original TITAN
#'   funtion call.
#' @param minSplt The minimum bin size for partitioning along the
#'   environmental gradient.  The default is to pass on the argument
#'   from the original TITAN funtion call.
#' @param nPerm The number of replicates used by the permutation
#'   procedure (not to be confused with the number of bootstrap
#'   replicates). The default is to pass on the argument from the
#'   original TITAN funtion call.
#' @param memory A logical indicating whether scratch files should
#'   be used to store temporary data in order to preserve RAM during
#'   bootstrapping of large data sets.  The default is to pass on
#'   the argument from the original TITAN funtion call.
#' @param imax A logical indicating whether taxon-specific change
#'   points should be determined by IndVal maxima or z-score maxima
#'   (as in Baker and King 2010). The default is to pass on the
#'   argument from the original TITAN funtion call.
#' @return A list of four elements: \itemize{ \item{bt.metrics }{A
#'   matrix with nrow equal to number of taxa where the first column
#'   is the bootstrapped IndVal or z score maximum, the second is
#'   the environmental value, the third is the indicator direction,
#'   and the fourth is the p value at that point.} \item{ivzs }{Z
#'   scores for all taxa across candidate change points in the
#'   replicate sample} \item{bsrti }{A sorted version of the
#'   bootstrapped environmental gradient} \item{rspdr }{Response
#'   direction (1 or 2) for all taxa across candidate change points
#'   in the replicate sample} }
#' @references Baker, ME and RS King.  2010. A new method for
#'   detecting and interpreting biodiversity and ecological
#'   community thresholds. Methods in Ecology and Evolution 1(1):
#'   25:37.
#' @seealso \code{\link{b.getivz}}, \code{\link{ivzsums}}
#' @author M. Baker and R. King
#' @keywords TITAN bootstrap
tboot <- function(bSeq, env, taxa, ivTot = ivTot, minSplt = minSplt,
  nPerm = nPerm, memory = memory, imax = imax) {

  print(bSeq)
  ## Create empty arrays for storing sum(z-) and sum(z+) across
  ## envcls by replicate
  boot.metrics <- rep(NA, length(bSeq))
  numUnit = length(env)
  rnum = runif(1, 0, 100)
  max.btSumz <- rep(NA, 2)

  for (i in 1:length(bSeq)) {
    ## Resample data with replacement, repeat part 1
    ## preliabiminaries
    nuid <- sample(numUnit, replace = TRUE)
    bTxa <- taxa[nuid, ]
    numTxa <- ncol(bTxa)
    bEnv <- env[nuid]
    bRankEnv <- rank(bEnv, ties.method = "random")
    bSrti <- sort(bEnv)
    bSrti2 <- sort(bRankEnv)
    bEnvCls <- bSrti[(minSplt):(numUnit - minSplt)]
    nClass <- length(bEnvCls)
    bEclass <- matrix(NA, numUnit, nClass)
    for (c in 1:nClass) {
      bEclass[, c] <- ((bRankEnv > bSrti2[minSplt + (c -
        1)]) * 1) + 1
    }
    boot.env <- bEclass
    boot.taxa <- bTxa

    ## Create temporary matrix for storing z scores and group
    ## assignments
    ivzScores.bt <- b.getivz(boot.env, boot.taxa, ivTot, nPerm = nPerm,
      imax = imax, numClass = nClass)

    ## Replace NaNs with 0s
    for (k in (numTxa + 1):(numTxa * 3)) {
      for (j in 1:nClass) {
        if (is.nan(ivzScores.bt[k, j])) {
          ivzScores.bt[k, j] <- 0
        }
      }
    }


    bt.metrics <- matrix(NA, numTxa, 4)  #output matrix of metrics

    ## Obtain env level for IndVal maxima and group assignment for
    ## each replicate bt.metric 1:4 include maxgrp, env.cp (either
    ## imax or zmax),z.score,iv.p
    for (j in 1:numTxa) {
      if (imax == FALSE) {
        maxSplt = which.max(abs(ivzScores.bt[j + (numTxa),
          ]))
      } else {
        maxSplt = which.max(abs(ivzScores.bt[j + (numTxa *
          2), ]))
      }
      bt.metrics[j, 1] <- ivzScores.bt[j, maxSplt]
      bt.metrics[j, 2] <- (bSrti[minSplt + maxSplt - 1] +
        bSrti[minSplt + maxSplt])/2
      bt.metrics[j, 3] <- ivzScores.bt[j + numTxa, maxSplt]
      bt.metrics[j, 4] <- ivzScores.bt[(j + (numTxa * 3)),
        maxSplt]
    }

    ## Obtain z scores and response direction
    ivzs <- ivzScores.bt[(numTxa + 1):(numTxa * 2), ]
    rspdr <- ivzScores.bt[1:numTxa, ]
    if (memory) {
      if (file.exists("temp.dir")) {
        write.table(ivzs, paste("temp.dir/boot.z.", rnum,
          ".", i, ".txt", sep = ""), append = F, row.names = F,
          col.names = F)
        write.table(rspdr, paste("temp.dir/boot.rd.", rnum,
          ".", i, ".txt", sep = ""), append = F, row.names = F,
          col.names = F)
        write.table(bSrti, paste("temp.dir/bSrti.", rnum,
          ".", i, ".txt", sep = ""), append = F, row.names = F,
          col.names = F)
      } else {
        dir.create("temp.dir")
        write.table(ivzs, paste("temp.dir/boot.z.", rnum,
          ".", i, ".txt", sep = ""), append = F, row.names = F,
          col.names = F)
        write.table(rspdr, paste("temp.dir/boot.rd.", rnum,
          ".", i, ".txt", sep = ""), append = F, row.names = F,
          col.names = F)
        write.table(bSrti, paste("temp.dir/bSrti.", rnum,
          ".", i, ".txt", sep = ""), append = F, row.names = F,
          col.names = F)
      }
      bt.list <- list(bt.metrics, rnum)
    } else {
      bt.list <- list(bt.metrics, ivzs, bSrti, rspdr)
    }
    boot.metrics[i] <- list(bt.list)
  }
  return(boot.metrics)
}






























#' Controls the allocation of bootstrap replicates
#'
#' A wrapper function for controlling the implementation of
#' bootstrap replicates using the function 'tboot' by sequential or
#' parallel processing.
#'
#' If 'ncpus'>1 evaluates to TRUE, the function employs the package
#' 'snow' to implement parallel processing on multicore processors
#' common in modern desktop computers.  With some minor modification
#' it is possible to configure this code to allocate processes to
#' cores on a high-performance computing cluster (i.e., a
#' supercomputer).  If 'ncpus'>1 evaluates to FALSE, the function
#' uses 'lapply' to run 'tboot' in sequence 1:nBoot times.
#'
#' @param env A vector of values for each sampling location along
#'   the environmental gradient.
#' @param taxa A site-by-taxon matrix of taxa counts at each
#'   sampling location.
#' @param ivTot A logical indicating whether IndVal scores should be
#'   calculated using total relative abundance or the mean relative
#'   abundace originally proposed by Dufrene and Legendre (1997).
#'   The default is to pass on the argument from the original TITAN
#'   funtion call.
#' @param boot A logical specifying whether or not to implement
#'   TITAN's' boostrap procedure. The default is to use the value
#'   specified in the original TITAN function call.
#' @param ncpus An argument specifying the number of processing
#'   cores used by the TITAN function call.  If ncpus>1 then
#'   parallel processing is implemented.  The default is to use the
#'   value specified in the original TITAN function call.
#' @param nBoot An argument specifying the number of bootstrap
#'   replicates.  The default is to use the value specified in the
#'   original TITAN function call.
#' @param minSplt An argument specifying minimum split size for
#'   partitioning along the environmental gradient.  The default is
#'   to use the value specified in the original TITAN function call.
#' @param nPerm The number of replicates used by the permutation
#'   procedure (not to be confused with the number of bootstrap
#'   replicates).
#' @param memory A logical indicating whether scratch files should
#'   be used to store temporary data in order to conserve active
#'   memory during bootstrapping of large data sets.  The default is
#'   to pass on the argument from the original TITAN funtion call.
#' @param imax A logical indicating whether taxon-specific change
#'   points should be determined by IndVal maxima or z-score maxima
#'   (as in TITAN v1.0). The default is to pass on the argument from
#'   the original TITAN funtion call.
#' @param numUnit An argument specifying the number of values along
#'   the environmental gradient.
#' @return A list of two items: \itemize{ \item{bSeq }{An index of
#'   the sequence of bootstrap replicates.  The structure of bSeq
#'   will differ for sequential or parallel processing.}
#'   \item{ivz.bt.list }{Itself a list of four items comprising
#'   output passed on from function 'tboot'} }
#' @references Baker, ME and RS King.  2010. A new method for
#'   detecting and interpreting biodiversity and ecological
#'   community thresholds. Methods in Ecology and Evolution 1(1):
#'   25:37.
#' @seealso \code{\link{tboot}}, \code{\link{small.boot}},
#'   \code{\link{big.boot}}, \code{\link{titan}}
#' @author M. Baker and R. King
#' @keywords TITAN bootstrap
boot.titan <- function(env, taxa, ivTot = ivTot, boot = boot, ncpus = ncpus,
  nBoot = nBoot, minSplt = minSplt, nPerm = 250, memory = memory,
  imax = imax, numUnit = numUnit) {

  ## If multiple cores are available, take advantage of parallel
  ## processing
  if (ncpus > 1) {
    print(paste("Bootstrap resampling in parallel using ",
      ncpus, " CPUs", "...no index will be printed to screen",
      sep = ""))
    cores <- rep("localhost", ncpus)
    cl <- parallel::makeCluster(cores, type = "SOCK")
    # parallel::clusterExport(cl,
    # c('indval','indvals','permiv','b.getivz','indvalp','indvalps','ivzsums'))
    bSeq = parallel::clusterSplit(cl, 1:nBoot)
    ivz.bt.list = parallel::clusterApply(cl, bSeq, tboot, env = env,
      taxa = taxa, ivTot = ivTot, minSplt = minSplt, nPerm = nPerm,
      memory = memory, imax = imax)
    parallel::stopCluster(cl)
  } else {
    ## otherwise run bootstrap in sequence
    print("Bootstrap resampling in sequence...")
    bSeq = 0
    ivz.bt.list = lapply(1:nBoot, tboot, env = env, taxa = taxa,
      ivTot = ivTot, minSplt = minSplt, nPerm = nPerm, memory = memory,
      imax = imax)
  }

  list(bSeq, ivz.bt.list)
}























#' Summarizes raw output from TITAN's bootstrap procedure
#'
#' A function to take output from TITAN's bootstrap procedure and
#' process it for summary output.  The default is to perform this
#' processing entirely within active memory, but in the event of
#' overflowing system capacity, an optional program writes temporary
#' files to a scratch directory to circumvent memory limits.
#'
#' Use of 'small.boot' versus 'big.boot' is controlled by the
#' argument 'memory' in the original TITAN function call and passed
#' to the wrapper function 'titan'.  The two progams have identical
#' functionality, but they accomplish those functions differently to
#' deal with memory limitations.
#'
#' For sequential processing of the bootsrtap, the index 'bSeq' is
#' simply a sequence from 1:nBoot that is printed to the screen. For
#' parallel processing, 'bSeq' is a list of length equal to 'ncpus',
#' where each item is a segment of the sequence allocated to each
#' processing core.  Thus, depending on whether 'ncpus'>1, the value
#' of 'bSeq' is used differently to extract values from the
#' bootstrap output list.
#'
#' The first part of each function consists of defining output
#' matrices, the second involves extraction of output from the
#' bootstrap list, the third part involves calculating purity,
#' reliability, the median z score, and quantiles of the
#' bootstrapped change points for each taxon.  These values are used
#' to complete the 'sppmax' output table and to identify the taxa
#' that meet purity and reliability criteria.  The final portion of
#' each function finds the maximum sum(z-), sum(z+), f.sum(z-), and
#' f.sum(z+) for each bootstrap replicate for later estimation of
#' confidence intervals.  The final portion of the summary involves
#' calculating the filtered and unfiltered sum(z) scores for each
#' bootstrap replicate from the matrix of z scores and response
#' directions passed from the function boot.titan() within
#' ivz.bt.list
#'
#' @param ivz.bt.list A list of output from each bootstrap replicate
#'   passed from \code{\link{boot.titan}}.
#' @param bSeq An index of the sequence of bootstrap replicates.
#' @param sppmax A taxon-specific summary output table for TITAN.
#' @param obs1 A binary vector indicating membership in the
#'   decreasing group of taxa.
#' @param obs2 A binary vector indicating membership in the
#'   increasing group of taxa.
#' @param nBoot An argument specifying the number of bootstrap
#'   replicates.  The default is to use the value specified in the
#'   original TITAN function call.
#' @param numClass An argument specifying the number of candidate
#'   partitions along the environmental gradient.
#' @param numUnit An argument specifying the number of values along
#'   the environmental gradient.
#' @param ncpus An argument specifying the number of processing
#'   cores used by the TITAN function call.  If ncpus>1 then
#'   parallel processing is indicated.  The default is to use the
#'   value specified in the original TITAN function call.
#' @param pur.cut An argument specifying the cutoff value for
#'   determining purity.  The default is to use the value specified
#'   in the original TITAN function call.
#' @param rel.cut An argument specifying the cutoff value for
#'   determining reliability.  The default is to use the value
#'   specified in the original TITAN function call.
#' @param minSplt An argument specifying minimum split size of
#'   partitioning along the environmental gradient.  The default is
#'   to use the value specified in the original TITAN function call.
#' @return A list of six items: \itemize{ \item{sppSub1 }{A vector
#'   of taxon index numbers for pure and reliable decreasers}
#'   \item{sppSub2 }{A vector of taxon index numbers for pure and
#'   reliable increasers} \item{sppmax }{The completed
#'   taxon-specific summary output table for TITAN} \item{maxSumz
#'   }{A 2-column matrix of environmental values at sum(z-) and
#'   sum(z+) maxima across all bootstrap replicates} \item{maxFsumz
#'   }{A 2-column matrix of environmental values at filtered sum(z-)
#'   and sum(z+) maxima across all bootstrap replicates}
#'   \item{metricArray}{An array of group membership, env change
#'   points, z scores, and p values for passing to 'plot.IVecdf'} }
#' @references Baker, ME and RS King.  2010. A new method for
#'   detecting and interpreting biodiversity and ecological
#'   community thresholds. Methods in Ecology and Evolution 1(1):
#'   25:37.
#' @references Baker ME and RS King. 2013. Of TITAN and straw men:
#'   an appeal for greater understanding of community data.
#'   Freshwater Science 32(2):489-506.
#' @seealso \code{\link{boot.titan}}, \code{\link{tboot}},
#'   \code{\link{titan}}
#' @author M. Baker and R. King
#' @keywords TITAN purity reliability sum(z)
#' @name smallBigBoot






#' @rdname smallBigBoot
small.boot <- function(ivz.bt.list, bSeq, sppmax, obs1, obs2, nBoot,
  numClass, numUnit, ncpus, pur.cut, rel.cut, minSplt) {

  numTxa <- nrow(sppmax)
  ## Create empty matrix for storing IndVal maxima and group
  ## assignments
  metricArray <- array(NA, c(numTxa, 4, nBoot))

  ## Create empty matrices for storing nuids and z values at
  ## IndVal maxima for bootreps
  zArray <- array(NA, c(numTxa, numClass, nBoot))
  bEnvMatrix <- matrix(NA, nBoot, numUnit)
  rspdir <- array(NA, c(numTxa, numClass, nBoot))

  ## Store env values of filtered and unfiltered sum(z) maxima by
  ## bootstrap replicate
  sumzBoot <- matrix(NA, nBoot, 2)
  sumzBoot.f <- matrix(NA, nBoot, 2)
  maxSumz <- matrix(NA, nBoot, 2)
  maxFsumz <- matrix(NA, nBoot, 2)
  aseq = 0

  ### The following loop creates an array across all bootreps from
  ### parallel or sequenced result list
  if (ncpus > 1) {
    for (s in 1:ncpus) {
      for (l in 1:length(bSeq[[s]])) {
        aseq = aseq + 1
        metricArray[, , aseq] <- unlist(ivz.bt.list[[s]][[l]][[1]])
        zArray[, , aseq] <- unlist(ivz.bt.list[[s]][[l]][[2]])
        bEnvMatrix[aseq, ] <- unlist(ivz.bt.list[[s]][[l]][[3]])
        rspdir[, , aseq] <- unlist(ivz.bt.list[[s]][[l]][[4]])
      }
    }
  } else {
    for (i in 1:nBoot) {
      metricArray[, , i] <- unlist(ivz.bt.list[[i]][[1]][[1]])
      zArray[, , i] <- unlist(ivz.bt.list[[i]][[1]][[2]])
      bEnvMatrix[i, ] <- unlist(ivz.bt.list[[i]][[1]][[3]])
      rspdir[, , i] <- unlist(ivz.bt.list[[i]][[1]][[4]])
    }
  }

  purity1 <- rowMeans(metricArray[, 1, ] == 1, na.rm = T)
  purity2 <- rowMeans(metricArray[, 1, ] == 2, na.rm = T)
  reliab <- rowMeans(metricArray[, 4, ] < 0.05, na.rm = T)
  z.median <- apply(metricArray[, 3, ], 1, median, na.rm = T)
  quantiles <- t(apply(metricArray[, 2, ], 1, quantile, probs = c(0.05,
    0.1, 0.5, 0.9, 0.95), na.rm = T))
  sppmax[, 14] <- reliab
  sppmax[, 15] <- z.median
  sppmax[which(purity1 >= pur.cut & reliab >= rel.cut & obs1),
    16] <- 1
  sppmax[which(purity2 >= pur.cut & reliab >= rel.cut & obs2),
    16] <- sppmax[which(purity2 >= pur.cut & reliab >= rel.cut &
    obs2), 16] + 2
  sppmax[which(sppmax[, 4] == 1), 13] <- purity1[which(sppmax[,
    4] == 1)]
  sppmax[which(sppmax[, 4] == 2), 13] <- purity2[which(sppmax[,
    4] == 2)]
  sppmax[, 8:12] <- quantiles
  sppSub1 <- which(sppmax[, 16] == 1)
  sppSub2 <- which(sppmax[, 16] == 2)

  for (i in 1:nBoot) {
    sumzBoot[i, 1] <- which.max(colSums(zArray[, , i] * (rspdir[,
      , i] == 1), na.rm = T))
    if (length(sppSub1) > 1) {
      sumzBoot.f[i, 1] <- which.max(colSums(zArray[sppSub1,
        , i] * (rspdir[sppSub1, , i] == 1), na.rm = T))
    } else {
      if (length(sppSub1) > 0) {
        sumzBoot.f[i, 1] <- which.max(zArray[sppSub1, ,
          i] * (rspdir[sppSub1, , i] == 1))
      }
    }
    sumzBoot[i, 2] <- which.max(colSums(zArray[, , i] * (rspdir[,
      , i] == 2), na.rm = T))
    if (length(sppSub2) > 1) {
      sumzBoot.f[i, 2] <- which.max(colSums(zArray[sppSub2,
        , i] * (rspdir[sppSub2, , i] == 2), na.rm = T))
    } else {
      if (length(sppSub2) > 0) {
        sumzBoot.f[i, 2] <- which.max(zArray[sppSub2, ,
          i] * (rspdir[sppSub2, , i] == 2))
      }
    }
    maxSumz[i, 1] <- (bEnvMatrix[i, (minSplt + sumzBoot[i,
      1]) - 1] + bEnvMatrix[i, (minSplt + sumzBoot[i, 1])])/2
    maxSumz[i, 2] <- (bEnvMatrix[i, (minSplt + sumzBoot[i,
      2]) - 1] + bEnvMatrix[i, (minSplt + sumzBoot[i, 2])])/2
    maxFsumz[i, 1] <- (bEnvMatrix[i, (minSplt + sumzBoot.f[i,
      1]) - 1] + bEnvMatrix[i, (minSplt + sumzBoot.f[i, 1])])/2
    maxFsumz[i, 2] <- (bEnvMatrix[i, (minSplt + sumzBoot.f[i,
      2]) - 1] + bEnvMatrix[i, (minSplt + sumzBoot.f[i, 2])])/2
  }
  list(sppSub1, sppSub2, sppmax, maxSumz, maxFsumz, metricArray)
}















#' @rdname smallBigBoot
big.boot <- function(ivz.bt.list, bSeq, sppmax, obs1, obs2, nBoot,
  numClass, numUnit, ncpus, pur.cut, rel.cut, minSplt) {

  numTxa <- nrow(sppmax)
  ## Create empty matrix for storing IndVal maxima and group
  ## assignments
  metricArray <- array(NA, c(numTxa, 4, nBoot))

  ## Store env values of filtered and unfiltered sum(z) maxima by
  ## bootstrap replicate
  sumzBoot <- matrix(NA, nBoot, 2)
  sumzBoot.f <- matrix(NA, nBoot, 2)
  maxSumz <- matrix(NA, nBoot, 2)
  maxFsumz <- matrix(NA, nBoot, 2)
  track.seq <- rep(NA, nBoot)
  aseq = 0

  ### The following loop creates an array across all bootreps from
  ### parallel or sequenced result list
  if (ncpus > 1) {
    for (s in 1:ncpus) {
      for (l in 1:length(bSeq[[s]])) {
        aseq = aseq + 1
        metricArray[, , aseq] <- unlist(ivz.bt.list[[s]][[l]][[1]])
        track.seq[aseq] <- ivz.bt.list[[s]][[l]][[2]]
      }
    }
  } else {
    for (i in 1:nBoot) {
      aseq = aseq + 1
      metricArray[, , i] <- unlist(ivz.bt.list[[i]][[1]][[1]])
      track.seq[aseq] <- ivz.bt.list[[i]][[1]][[2]]
    }
  }

  purity1 <- rowMeans(metricArray[, 1, ] == 1, na.rm = T)
  purity2 <- rowMeans(metricArray[, 1, ] == 2, na.rm = T)
  reliab <- rowMeans(metricArray[, 4, ] < 0.05, na.rm = T)
  z.median <- apply(metricArray[, 3, ], 1, median, na.rm = T)
  quantiles <- t(apply(metricArray[, 2, ], 1, quantile, probs = c(0.05,
    0.1, 0.5, 0.9, 0.95), na.rm = T))

  sppmax[, 14] <- reliab
  sppmax[, 15] <- z.median
  sppmax[which(purity1 >= pur.cut & reliab >= rel.cut & obs1),
    16] <- 1
  sppmax[which(purity2 >= pur.cut & reliab >= rel.cut & obs2),
    16] <- sppmax[which(purity2 >= pur.cut & reliab >= rel.cut &
    obs2), 16] + 2
  sppmax[which(sppmax[, 4] == 1), 13] <- purity1[which(sppmax[,
    4] == 1)]
  sppmax[which(sppmax[, 4] == 2), 13] <- purity2[which(sppmax[,
    4] == 2)]
  sppmax[, 8:12] <- quantiles
  sppSub1 <- which(sppmax[, 16] == 1)
  sppSub2 <- which(sppmax[, 16] == 2)

  ## Read z score matrices in one at a time and obtain filtered
  ## sumz maxima
  if (ncpus > 1) {
    aseq = 0
    for (s in 1:ncpus) {
      for (l in 1:length(bSeq[[s]])) {
        aseq = aseq + 1
        ztab = read.table(paste("temp.dir/boot.z.", track.seq[aseq],
          ".", l, ".txt", sep = ""))
        rspdir = read.table(paste("temp.dir/boot.rd.",
          track.seq[aseq], ".", l, ".txt", sep = ""))
        bSrtEnv = read.table(paste("temp.dir/bSrti.", track.seq[aseq],
          ".", l, ".txt", sep = ""))
        sumzBoot[aseq, 1] <- which.max(colSums(ztab * (rspdir ==
          1), na.rm = T))

        if (length(sppSub1) > 1) {
          sumzBoot.f[aseq, 1] <- which.max(colSums(ztab[sppSub1,
          ] * (rspdir[sppSub1, ] == 1), na.rm = T))
        } else {
          sumzBoot.f[aseq, 1] <- which.max(ztab[sppSub1,
          ] * (rspdir[sppSub1, ] == 1))
        }
        sumzBoot[aseq, 2] <- which.max(colSums(ztab * (rspdir ==
          2), na.rm = T))
        if (length(sppSub2) > 1) {
          sumzBoot.f[aseq, 2] <- which.max(colSums(ztab[sppSub2,
          ] * (rspdir[sppSub2, ] == 2), na.rm = T))
        } else {
          sumzBoot.f[aseq, 2] <- which.max(ztab[sppSub2,
          ] * (rspdir[sppSub2, ] == 2))
        }

        maxSumz[aseq, 1] <- (bSrtEnv[(minSplt + sumzBoot[aseq,
          1]) - 1, ] + bSrtEnv[(minSplt + sumzBoot[aseq,
          1]), ])/2
        maxSumz[aseq, 2] <- (bSrtEnv[(minSplt + sumzBoot[aseq,
          2]) - 1, ] + bSrtEnv[(minSplt + sumzBoot[aseq,
          2]), ])/2
        maxFsumz[aseq, 1] <- (bSrtEnv[(minSplt + sumzBoot.f[aseq,
          1]) - 1, ] + bSrtEnv[(minSplt + sumzBoot.f[aseq,
          1]), ])/2
        maxFsumz[aseq, 2] <- (bSrtEnv[(minSplt + sumzBoot.f[aseq,
          2]) - 1, ] + bSrtEnv[(minSplt + sumzBoot.f[aseq,
          2]), ])/2
      }
    }
  } else {
    for (i in 1:nBoot) {
      ztab = read.table(paste("temp.dir/boot.z.", track.seq[i],
        ".", 1, ".txt", sep = ""))
      rspdir = read.table(paste("temp.dir/boot.rd.", track.seq[i],
        ".", 1, ".txt", sep = ""))
      bSrtEnv = read.table(paste("temp.dir/bSrti.", track.seq[i],
        ".", 1, ".txt", sep = ""))
      sumzBoot[i, 1] <- which.max(colSums(ztab * (rspdir ==
        1), na.rm = T))
      sumzBoot[i, 2] <- which.max(colSums(ztab * (rspdir ==
        2), na.rm = T))
      sumzBoot.f[i, 1] <- which.max(colSums(ztab[sppSub1,
        , i] * (rspdir[sppSub1, , i] == 1), na.rm = T))
      sumzBoot.f[i, 2] <- which.max(colSums(ztab[sppSub2,
        , i] * (rspdir[sppSub2, , i] == 2), na.rm = T))
      maxSumz[i, 1] <- (bSrtEnv[(minSplt + sumzBoot[i, 1]) -
        1, ] + bSrtEnv[(minSplt + sumzBoot[i, 1]), ])/2
      maxSumz[i, 2] <- (bSrtEnv[(minSplt + sumzBoot[i, 2]) -
        1, ] + bSrtEnv[(minSplt + sumzBoot[i, 2]), ])/2
      maxFsumz[i, 1] <- (bSrtEnv[(minSplt + sumzBoot.f[i,
        1]) - 1, ] + bSrtEnv[(minSplt + sumzBoot.f[i, 1]),
        ])/2
      maxFsumz[i, 2] <- (bSrtEnv[(minSplt + sumzBoot.f[i,
        2]) - 1, ] + bSrtEnv[(minSplt + sumzBoot.f[i, 2]),
        ])/2
    }
  }
  list(sppSub1, sppSub2, sppmax, maxSumz, maxFsumz, metricArray)
}
