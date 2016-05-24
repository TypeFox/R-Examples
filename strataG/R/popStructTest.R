#' @name popStructTest 
#' @title Population Differentiation Tests
#' @description Conduct overall and/or pairwise tests of 
#'   population differentiation.
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param nrep number specifying number of permutation replicates to use for 
#'   permutation test.
#' @param stats a character vector or list of functions specifying which 
#'   anlayses to conduct. If characters, then valid possible choices are: 
#'   "phist", "fst", "fst.prime", "fis", "gst", "gst.prime", "gst.dbl.prime", "d", 
#'   or "chi2", or "all". If a list, then functions must be a valid population 
#'   structure function (see \code{\link{popStructStat}}) taking a 
#'   \linkS4class{gtypes} object and returning a named statistic estimate.
#' @param type character determining type of test to conduct. Can be "overall", 
#'   "pairwise", or "both". If "pairwise" or "both" are chosen and there are 
#'   only two strata, then only an overall test will be conducted.
#' @param max.cores The maximum number of cores to use to distribute separate 
#'   statistics over. Default (NULL) sets value to what is reported by 
#'   \code{\link[parallel]{detectCores} - 1}. Any value greater than this will 
#'   be set to this value. If \code{detectCores} reports \code{NA}, 
#'   \code{max.cores} will be set to 1.
#' @param keep.null logical. Keep the null distribution from the 
#'   permutation test?
#' @param quietly logical. Print progress and results?
#' @param write.output logical. Write a .csv file with results?
#' @param ... other parameters to be passed to population 
#'   differentiation functions.
#' 
#' @return
#' \describe{
#'  \item{overall}{a list containing:
#'    \tabular{ll}{
#'      \code{strata.freq} \tab a vector of the sample sizes for each stratum.\cr
#'      \code{result} \tab a matrix with the statistic estimate and p-value 
#'        for each statistic.\cr
#'      \code{null.dist} \tab a matrix with the null distributions for 
#'        each statistic.\cr
#'    }}
#'  \item{pairwise}{a list containing:
#'    \tabular{ll}{
#'      \code{result} \tab a data.frame with the result of each pairwise 
#'        comparison on a line.\cr
#'      \code{pair.mat} \tab a list with a pairwise matrix for each statistic. 
#'        Values in lower left are the statistic estimate, and upper right are p-values.\cr
#'      \code{null.dist} \tab a matrix with the null distributions for 
#'        each statistic.\cr
#'    }}
#' }
#' 
#' @note On multi-core systems, runs of separate statistics are automatically 
#'   distributed over as many cores as available (minus one). This can be 
#'   controlled by the \code{max.cores} argument if less core usage is 
#'   desired. 
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples 
#' data(msats.g)
#' msats.g <- stratify(msats.g, "fine")
#' 
#' # Just an overall Chi-squared test
#' ovl <- overallTest(msats.g, stats = "chi2", nrep = 100)
#' ovl
#' 
#' #' Just a pairwise test for Gst
#' pws <- pairwiseTest(msats.g, stats = list(statGst), nrep = 100)
#' pws
#' 
#' \dontrun{
#' #' Both overall and pairwise tests for Fst and F'st
#' full <- popStructTest(msats.g, stats = c("fst", "fst.prime"))
#' print(full$overall)
#' print(full$pairwise)
#' }
#' 
#' @importFrom utils write.csv
#' @export
#' 
popStructTest <- function(g, nrep = 1000, stats = "all", 
                          type = c("both", "overall", "pairwise"),
                          keep.null = FALSE, quietly = FALSE, 
                          max.cores = NULL, write.output = FALSE, ...) {
  # check arguments
  type <- match.arg(type)
    
  # conduct overall test
  overall <- NULL
  if(type %in% c("both", "overall")) {
    overall <- overallTest(
      g = g, nrep = nrep, stats = stats, keep.null = keep.null, 
      quietly = quietly, max.cores = max.cores, ...
    )    
  }
  
  # conduct pairwise test
  pairwise <- NULL
  if(type %in% c("both", "pairwise") & nStrata(g) > 2) {
    pairwise <- pairwiseTest(
      g = g, nrep = nrep, stats = stats, keep.null = keep.null, 
      quietly = quietly, max.cores = max.cores, ...
    )
  }

  if(write.output) {
    if(!is.null(overall)) {
      out.file <- gsub("[[:punct:]]", ".", 
       paste(description(g), "permutation test results.csv")
      )
      write.csv(overall$result, out.file)
    }
    if(!is.null(pairwise)) {
      for(stat in names(pairwise$pair.mat)) {
        out.file <- gsub("[[:punct:]]", ".", 
          paste(description(g), stat, "pairwise matrix.csv")
        )
        write.csv(pairwise$pair.mat[[stat]], out.file)
      }
    }    
  }
    
  invisible(list(overall = overall, pairwise = pairwise))
}


#' @rdname popStructTest
#' @importFrom parallel parLapply stopCluster detectCores
#' @importFrom swfscMisc pVal
#' @export
#' 
overallTest <- function(g, nrep = 1000, stats = "all", keep.null = FALSE, 
                        quietly = FALSE, max.cores = NULL, ...) {  
  
  stat.list <- statList(stats)
  if(length(stat.list) == 0) stop("no stats specified. NULL returned.")
  if(is.null(max.cores)) max.cores <- detectCores() - 1
  if(is.na(max.cores)) max.cores <- 1
  if(max.cores < 1) max.cores <- 1
  num.cores <- min(length(stat.list), max.cores)
  
  # check replicates
  if(is.null(nrep)) nrep <- 0
  if(!is.numeric(nrep) & length(nrep) != 1) {
    stop("'nrep' must be a single-element numeric vector")
  }
  if(nrep < 1) keep.null <- FALSE
  
  # remove unstratified samples
  if(any(is.na(strata(g)))) g <- g[, , strataNames(g)]
  
  # delete loci with no genotypes in at least one stratum
  to.delete <- unique(unlist(lapply(strataSplit(g), function(st.g) {
    n.genotyped <- nInd(st.g) - numMissing(st.g)
    names(n.genotyped)[n.genotyped == 0]
  })))
  if(length(to.delete) > 0) {
    warning(paste(
      "The following loci will be removed because they have no genotypes in one or more strata: ",
      paste(to.delete, collapse = ", ")
    ))
    g <- g[, setdiff(locNames(g), to.delete), ]
  }
  
  if(nStrata(g) == 1) stop("'g' must have more than one stratum defined.")
  
  if(!quietly) cat(
    cat("\n<<<", description(g), ">>>\n"),
    format(Sys.time()), ": Overall test :", nrep, "permutations\n"
  )
  
  # run each statistic and store results in a list
  strata.mat <- .permStrata(g, nrep) 
  stat.func <- function(f, strata.mat, keep.null, ...) {
    f(g, strata.mat = strata.mat, keep.null = keep.null, ...)
  }
  result <- if(num.cores == 1) {
    lapply(stat.list, stat.func, strata.mat = strata.mat, keep.null = keep.null, ...)
  } else {
    cl <- .setupClusters(num.cores)
    tryCatch({
      parLapply(cl, stat.list, stat.func, strata.mat = strata.mat, keep.null = keep.null, ...)
    }, finally = stopCluster(cl))
  }
  
  # create matrix of estimates and p-values
  result.mat <- t(sapply(result, function(x) x$result))
  rownames(result.mat) <- sapply(result, function(x) x$stat.name)
  
  # collect null distribution
  null.dist <- if(keep.null) {
    nd <- sapply(result, function(x) x$null.dist)
    colnames(nd) <- rownames(result.mat)
    nd
  } else NULL
  
# # calculate list of observed values for each population structure function
#   # conduct permutation test
#   null.dist <- NULL
#   if(nrep > 0 & length(stat.list) > 0) {
#     st <- lapply(1:nrep, function(i) sample(strata(g)))  
#     
#     # setup permutation function to return vector of values from each population 
#     #   structure statistic in stat.funcs
#     perm.func <- function(ran.strata, g, stat.funcs, ...) {
#       sapply(1:length(stat.funcs), function(j) {
#         stat.funcs[[j]](g, strata = ran.strata, ...)
#       })
#     }
#     
#     # collect null distribution
#     if(num.cores > 1) {
#       # setup clusters
#       cl <- .setupClusters(num.cores)
#       tryCatch({
#         # calculate matrix of null distributions
#         null.dist <- parLapply(cl, st, perm.func, g = g, stat.funcs = stat.list, ...)
#       })
#       stopCluster(cl)
#       closeAllConnections()
#     } else {
#       null.dist <- lapply(st, perm.func, g = g, stat.funcs = stat.list, ...)
#     }
#     null.dist <- do.call(rbind, null.dist)
#     colnames(null.dist) <- rownames(result)
#     
#     # calculate vector of p-values
#     for(x in rownames(result)) {
#       est <- result[x, "estimate"]
#       if(!is.na(est)) result[x, "p.val"] <- pVal(est, na.omit(null.dist[, x]))
#     }
#   } 
  # if(!keep.null) null.dist <- NULL
  
  # collect strata frequencies to named vector 
  strata.freq <- table(strata(g), useNA = "no")
  
  if(!quietly) {
    cat("\n")
    print(cbind(N = strata.freq))
    cat("\nPopulation structure results:\n")
    print(result.mat)
    cat("\n")
  }
  
  invisible(
    list(strata.freq = strata.freq, result = result.mat, null.dist = null.dist)
  )
}


#' @rdname popStructTest
#' @export
#' 
pairwiseTest <- function(g, nrep = 1000, stats = "all", keep.null = FALSE, 
                         quietly = FALSE, max.cores = NULL, ...) { 
  
  if(nStrata(g) == 1) stop("'g' must have more than one stratum defined.")
  
  if(!quietly) cat(
    cat("\n<<<", description(g), ">>>\n"),
    format(Sys.time()), ": Pairwise tests :", nrep, "permutations\n"
  )
  
  # create strata pairs
  strata.pairs <- .strataPairs(g)
  
  # run permutation test on all pairwise gtypes subsets
  pair.list <- vector("list", length = nrow(strata.pairs))
  for(i in 1:nrow(strata.pairs)) {
    pair <- unlist(strata.pairs[i, ])
    
    if(!quietly) {
      cat("  ", format(Sys.time()), ":", paste(pair, collapse = " v. "), "\n")
    }
    
    pair.list[[i]] <- overallTest(
      g = g[ , , pair], nrep = nrep, stats = stats, keep.null = keep.null, 
      quietly = TRUE, max.cores = max.cores, ...
    )
  }
  
  # compile results in 'pair.list' into a data.frame
  result <- do.call(rbind, lapply(pair.list, function(pair) {
    result.vec <- as.vector(t(pair$result))
    names(result.vec) <- paste(
      rep(rownames(pair$result), each = 2), c("", ".p.val"), sep = ""
    ) 
    result.vec <- rbind(result.vec)
    s1 <- names(pair$strata.freq)[1]
    s2 <- names(pair$strata.freq)[2]
    n1 <- pair$strata.freq[1]
    n2 <- pair$strata.freq[2]
    strata.1 <- paste(s1, " (", n1, ")", sep = "")
    strata.2 <- paste(s2, " (", n2, ")", sep = "")
    df <- data.frame(pair.label = paste(strata.1, " v. ", strata.2, sep = ""), 
      strata.1 = s1, strata.2 = s2, n.1 = n1, n.2 = n2,
      stringsAsFactors = FALSE
    )
    cbind(df, result.vec)   
  }))
  rownames(result) <- NULL
  
  # create pairwise matrices - lower left is estimate, upper right is p-value 
  stat.cols <- seq(6, ncol(result), 2)
  strata <- sort(levels(strata(g)))
  mat <- matrix(nrow = length(strata), ncol = length(strata), 
    dimnames = list(strata, strata)
  )
  pair.mat <- lapply(stat.cols, function(i) {
    for(j in 1:nrow(result)) {
      strata.1 <- result$strata.1[j]
      strata.2 <- result$strata.2[j]
      mat[strata.2, strata.1] <- result[j, i]
      mat[strata.1, strata.2] <- result[j, i + 1]
    }
    mat
  })
  names(pair.mat) <- colnames(result)[stat.cols]
  
  # compile null distributions into list of matrices
  null.dist <- if(keep.null) {
    null.mat <- lapply(pair.list, function(pair) pair$null.dist)
    names(null.mat) <- result$pair.label
    null.mat
  } else NULL
  
  if(!quietly) {
    cat("\nPopulation structure results:\n")
    print(result[, c(1, 6:ncol(result))])
    cat("\n")
  }
  
  invisible(list(result = result, pair.mat = pair.mat, null.dist = null.dist))
}


#' @rdname popStructTest
#' @export
#' 
statList <- function(stats = "all") {
  # check stats and return list of functions
  stat.list <- list(
    chi2 = statChi2,
    d = statJostD,
    fst = statFst,
    fst.prime = statFstPrime,
    fis = statFis,
    gst = statGst,
    gst.prime = statGstPrime,
    gst.dbl.prime = statGstDblPrime,
    phist = statPhist
  )
  
  if(is.character(stats)) {
    stats <- tolower(stats)
    if("all" %in% stats) {
      stat.list 
    } else {
      missing <- !stats %in% names(stat.list)
      if(sum(missing) > 0) {
        missing <- paste(stats[missing], collapse = ", ")
        stop(paste("the following stats could not be found:", missing))
      }
      stat.list[stats]
    } 
  } else if(is.list(stats) & all(sapply(stats, is.function))) {
    stats
  } else {
    stop("'stats' is not a list of functions.")
  }
}