#' The primary function call for package TITAN
#'
#' The primary wrapper function controlling operation of all
#' subfunctions (txa.screen, env.part, getivz, ivzsums, obs.summ,
#' boot.titan, small/big.boot, sumz.tab) apart from plotting
#' functions within TITAN.
#'
#' @param env A vector of environmental values.
#' @param txa A site by taxon matrix containing observed counts at
#'   each sampling location.
#' @param minSplt The  minimum split size to be used in
#'   partitioning.
#' @param numPerm The number of replicates to be used during
#'   permutation.
#' @param boot A logical indicating whether bootstrap resampling
#'   should be performed.
#' @param nBoot The number of replicates to be used during bootstrap
#'   resampling.
#' @param imax A logical indication whether taxon-specific change
#'   points should be determined using IndVal maxima or z-score
#'   maxima (as in Baker and King 2010, 2013).
#' @param ivTot A logical indicating whether IndVal scores should be
#'   calculated using total relative abundance or the mean relative
#'   abundace originally proposed by Dufrene and Legendre (1997).
#' @param pur.cut A proportion specifying the cutoff value for
#'   determining purity across all bootstrap replicates.
#' @param rel.cut A proportion specifying the cutoff value for
#'   determining reliability across all bootstrap replicates.
#' @param ncpus The number of processing cores to be used during
#'   processing.  If greater than 1, TITAN will load and use the
#'   package 'snow' to perform parallel processing on each core.
#' @param memory A logical indicating whether temporary files should
#'   be written to a scrath directory during bootstrap processing to
#'   preserve active memory.  This fucntion is sometimes necessary
#'   for large data files (e.g., >400 sampling sites and >100 taxa).
#' @return A list with 13 items: \itemize{ \item{sppmax
#'   }{Description of 'comp1'} \item{sumz.cp }{Description of
#'   'comp1'} \item{env }{The vector of environmental values used in
#'   the TITAN function call} \item{taxa }{The site-by-taxon matrix
#'   used in the TITAN function call} \item{envlcs }{A vector of
#'   candidate partitions derived from subtracting 'minSplt' from
#'   'env'} \item{srtEnv }{A sorted version of environmental values}
#'   \item{ivzScores }{A matrix containing group membership, z
#'   scores, IndVals, and p values for each taxon at every candidate
#'   partition in 'envcls'} \item{ivz }{A 2-column matrix containing
#'   parallel vectors of sum(z-) and sum(z+ scores for every
#'   candidate partition in 'envcls')} \item{ivz.f }{A 2-column
#'   matrix containing parallel vectors of sum(z-) and sum(z+ scores
#'   filtered by pure and reliable taxa for every candidate
#'   partition in 'envcls')} \item{maxSumz }{A 2-column matrix of
#'   environmental values at sum(z-) and sum(z+) maxima across all
#'   bootstrap replicates} \item{maxFsumz }{A 2-column matrix of
#'   environmental values at filtered sum(z-) and sum(z+) maxima
#'   across all bootstrap replicates} \item{metricArray }{An array
#'   of group membership, env change points, z scores, and p values
#'   equivalent to 'ivzScores' for each bootstrap replicate}
#'   \item{targs }{A vector of arguments used in the TITAN function
#'   call} }
#' @references Baker, ME and RS King.  2010. A new method for
#'   detecting and interpreting biodiversity and ecological
#'   community thresholds.  Methods in Ecology and Evolution 1(1):
#'   25:37.
#' @references King, RS and ME Baker  2010. Considerations for
#'   identifying and interpreting ecological community thresholds.
#'   Journal of the North American Benthological Association
#'   29(3):998-1008.
#' @references Baker ME and RS King. 2013. Of TITAN and straw men:
#'   an appeal for greater understanding of community data.
#'   Freshwater Science 32(2):489-506.
#' @author M. Baker and R. King
#' @keywords TITAN sum(z) permutation bootstrap IndVal
#' @export
#' @examples
#'
#' \dontrun{ resampling is very slow
#'
#' data(glades.env)
#' data(glades.taxa)
#' glades.titan <- titan(glades.env, glades.taxa, minSplt = 5,
#'   numPerm = 250, boot = TRUE, nBoot = 500, imax = FALSE,
#'   ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95, ncpus = 1, memory = FALSE
#' )
#'
#' }
#'
titan <- function(env, txa, minSplt = 5, numPerm = 250, boot = TRUE,
  nBoot = 500, imax = FALSE, ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95,
  ncpus = 1, memory = FALSE) {

  # --------------------------------------------------------------------------
  # PART 1: Data Screen
  # --------------------------------------------------------------------------

  ## Load 'parallel' library for parallel processing
  ## if(ncpus>1){library(parallel)} Not necessary in package build
  ## because it is already imported.

  ## Run Preliminary Checks and Input Filters

  ## Abundance matrix (txa) should contain samples(rows) x
  ## species(columns)
  taxa <- txa.screen(txa, minSplt = minSplt)

  # if env is a data frame of one column, convert to vector
  if(is.data.frame(env) && ncol(env) == 1) env <- env[[1]]

  ## Environmental vector (env) should have same length as samples
  ## Sort sites by env, ensure candidate threshold group size
  ## always >= minsplit
  e.prt <- env.part(env, taxa, minSplt = minSplt)
  env <- e.prt[[1]]
  numUnit <- e.prt[[2]]
  numTxa <- e.prt[[3]]
  numClass <- e.prt[[4]]
  srtEnv <- e.prt[[5]]
  envcls <- e.prt[[6]]
  eclass <- e.prt[[7]]

  # --------------------------------------------------------------------------
  # PART 2: Obtain IndVals, z-scores, and sum(z) values
  # --------------------------------------------------------------------------

  ## Run IndVal calculation and permute for z score computations,
  ## track processing time
  ivzScores <- getivz(eclass, taxa, ivTot, numPerm, numClass)

  ## Sum IndVal z scores across all taxa
  ivz <- ivzsums(ivzScores)
  print("IndVal z score calculation complete")

  ## Summarize Observed Results
  obs.summary <- obs.summ(ivzScores, taxa, srtEnv, minSplt = minSplt,
    imax = imax)
  obs1 = obs.summary[[1]]
  obs2 = obs.summary[[2]]
  sppmax = obs.summary[[3]]
  rm(obs.summary)

  # --------------------------------------------------------------------------
  # PART 3: Bootstrap procedure for TITAN
  # --------------------------------------------------------------------------

  ## Run bootstrap procedure
  if (boot) {

    ## Control of Bootstrap Processing (Sequential or Parallel)
    ivzBtSeq <- boot.titan(env, taxa, ivTot, boot, ncpus, nBoot,
      minSplt, numPerm, memory, imax, numUnit)
    bSeq <- ivzBtSeq[[1]]
    ivz.bt.list <- ivzBtSeq[[2]]
    rm(ivzBtSeq)

    ## Summarize output, if file size too large write some to
    ## temporary storage on disk
    if (memory) {
      boot.summ <- big.boot(ivz.bt.list, bSeq, sppmax, obs1,
        obs2, nBoot, numClass, numUnit, ncpus, pur.cut,
        rel.cut, minSplt)
    } else {
      boot.summ <- small.boot(ivz.bt.list, bSeq, sppmax,
        obs1, obs2, nBoot, numClass, numUnit, ncpus, pur.cut,
        rel.cut, minSplt)
    }
    sppSub1 <- boot.summ[[1]]
    sppSub2 <- boot.summ[[2]]
    sppmax <- boot.summ[[3]]
    maxSumz <- boot.summ[[4]]
    maxFsumz <- boot.summ[[5]]
    metricArray <- boot.summ[[6]]
  } else {
    maxSumz <- 0
    ivz.f <- 0
    maxFsumz <- 0
    sppSub1 <- 0
    sppSub2 <- 0
    metricArray <- 0
    z.median <- 0
  }


  # --------------------------------------------------------------------------
  # PART 4: Summarize sumz TITAN output
  # --------------------------------------------------------------------------
  nsumz1 <- sum(sppmax[, 16] == 1)
  nsumz2 <- sum(sppmax[, 16] == 2)
  print(paste("Proportion of pure and reliable taxa = ", sum(nsumz2,
    nsumz1, na.rm = T)/nrow(sppmax), sep = ""))
  sumz.list <- sumz.tab(ivzScores, ivz, srtEnv, sppmax, maxSumz,
    maxFsumz, minSplt, boot)
  sumz.cp <- sumz.list[[1]]
  ivz.f <- sumz.list[[2]]

  if (nsumz1 > 2 & nsumz2 > 2) {
    print(sppmax)
    print(sumz.cp)
    print("TITAN complete")
  } else {
    print(sppmax)
    print("Warning: low number of pure and reliable taxa, sum(z) output should be interpreted with caution")
    print(paste("Number of z- taxa = ", nsumz1, ", Number of z+ taxa = ",
      nsumz2, sep = ""))
    print("TITAN complete")
  }

  # --------------------------------------------------------------------------
  # PART 5: Define output lists for TITAN object
  # --------------------------------------------------------------------------
  targs <- c(minSplt, numPerm, boot, nBoot, imax, ivTot, pur.cut,
    rel.cut, ncpus, memory)
  names(targs) <- c("minSplt", "numPerm", "boot", "nBoot", "imax",
    "ivTot", "pur.cut", "rel.cut", "ncpus", "memory")

  titan.out <- list(sppmax, sumz.cp, env, taxa, envcls, srtEnv,
    ivzScores, ivz, ivz.f, maxSumz, maxFsumz, metricArray,
    targs)
  names(titan.out) <- c("sppmax", "sumz.cp", "env", "taxa", "envcls",
    "srtEnv", "ivzScores", "ivz", "ivz.f", "maxSumz", "maxFsumz",
    "metricArray", "arguments")

  return(titan.out)


  ## END

}
