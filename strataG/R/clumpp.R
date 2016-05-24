#' @title Run CLUMPP
#' @description Run CLUMPP to aggregate multiple STRUCTURE runs.
#' 
#' @param sr result from \code{\link{structure}} or folder name containing 
#'   STRUCTURE output files.
#' @param k choice of \emph{k} in \code{sr} to combine.
#' @param align.algorithm algorithm to be used for aligning the runs. Can be 
#'   "full.search", "greedy", or "large.k".
#' @param sim.stat pairwise matrix similarity statistic to be used. Can be 
#'   "g" or "g.prime".
#' @param greedy.option input order of runs to be tested. Required if 
#'   \code{align.algorithm} is "greedy" or "large.k". Valid choices are:
#'   \tabular{ll}{
#'     \code{all} \tab test all possible input orders of runs (note that this 
#'       option increases the run-time sub-stantially unless R is small).\cr
#'     \code{ran.order} \tab test a specified number of random input orders of 
#'       runs set by the \code{repeats} parameter.\cr
#'   }
#' @param repeats the number of input orders of runs to be tested. Only used if 
#'   \code{align.algorithm} is "greedy" or "large.k", and \code{greedy.option} 
#'     is "ran.order".
#' @param order.by.run permute the clusters according to the cluster order of 
#'   a specific run. Set this parameter to a number from 1 to the number of 
#'   runs in \code{sr}.
#' @param label label to use for input and output files.
#' @param delete.files logical. Delete all files when CLUMPP is finished?
#' 
#' @note CLUMPP is not included with \code{strataG} and must be downloaded 
#'   separately. Additionally, it must be installed such that it can be run from 
#'   the command line in the current working directory. See the vignette 
#'   for \code{external.programs}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @references Mattias Jakobsson and Noah A. Rosenberg. 2007. CLUMPP: a cluster 
#'   matching and permutation program for dealing with label switching and 
#'   multimodality in analysis of population structure. 
#'   Bioinformatics 23(14):1801-1806. Available at 
#'   \url{http://web.stanford.edu/group/rosenberglab/clumppDownload.html}
#' 
#' @seealso \code{\link{structure}}
#' 
#' @importFrom utils file_test write.table
#' @export
#' 
clumpp <- function(sr, k, align.algorithm = "greedy", sim.stat = "g",
                      greedy.option = "ran.order", repeats = 100, 
                      order.by.run = 0, label = NULL, delete.files = TRUE) {
  
  if(k < 2) stop("k must be greater than 1.")
  
  if(!tolower(align.algorithm) %in% c("full.search", "greedy", "large.k")) {
    stop("'align.algorithm' must be either 'full.search', 'greedy', or 'large.k'.")
  }
  if(!tolower(sim.stat) %in% c("g", "g.prime")) {
    stop("'sim.stat' must be either 'g' or 'g.prime'.")
  }
    
  if(!"structure.result" %in% class(sr)) {
    folder <- sr
    if(!is.null(folder)) {
      if(!file.exists(folder)) dir.create(folder)
      if(!file_test("-d", folder)) stop("'folder' is not a valid folder.")
      label <- file.path(folder, label)
    }
  }
  
  if(!is.null(label)) label <- gsub(" ", ".", label)  
  param.file <- ifelse(is.null(label), "paramfile", 
                       paste(label, "paramfile", sep = "_"))
  ind.file <- ifelse(is.null(label), "indfile", 
                     paste(label, "indfile", sep = "_"))
  permutation.file <- ifelse(is.null(label), "permutationfile", 
                             paste(label, "permutationfile", sep = "_"))
  out.file <- ifelse(is.null(label), "outfile", 
                     paste(label, "outfile", sep = "_"))
  misc.file <- ifelse(is.null(label), "miscfile", 
                      paste(label, "miscfile", sep = "_"))
  
  eq.k <- sapply(sr, function(x) x$summary["k"] == k)
  if(sum(eq.k) == 0) stop(paste("no entries for k =", k, "found in 'sr'."))
  sr <- sr[eq.k]
  
  align.algorithm <- switch(align.algorithm, full.search = 1, 
                            greedy = 2, large.k = 3)
  sim.stat <- switch(sim.stat, g = 1, g.prime = 2)
  
  greedy.params <- if(align.algorithm != 1) {
    greedy.option <- switch(greedy.option, all = 1, ran.order = 2)
    gp <- paste("GREEDY_OPTION", greedy.option)
    if(greedy.option != 1) {
      if(is.null(repeats)) {
        stop("'repeats' must be specified if 'greedy.options' is 'greedy' or 'large.k'.")
      }
      gp <- c(gp, paste("REPEATS", repeats), 
              paste("PERMUTATIONFILE", permutation.file))
      # write permutationfile
      perm.mat <- t(sapply(1:repeats, function(i) sample(1:length(sr))))
      write.table(perm.mat, file = permutation.file, quote = FALSE, sep = " ", 
                  row.names = FALSE, col.names = FALSE)
    }
    gp
  } else NULL
  
  # write paramfile
  param.txt <- c(
    "DATATYPE 0",
    paste("INDFILE", ind.file),
    paste("OUTFILE", out.file),
    paste("MISCFILE", misc.file),
    paste("K", k),
    paste("C", nrow(sr[[1]]$q.mat)),
    paste("R", length(sr)),
    paste("M", align.algorithm),
    "W 0",
    paste("S", sim.stat),
    greedy.params,
    "PRINT_PERMUTED_DATA 0",
    "PRINT_EVERY_PERM 0",
    "PRINT_RANDOM_INPUTORDER 0",
    "OVERRIDE_WARNINGS 1",
    paste("ORDER_BY_RUN", order.by.run)
  )
  param.txt <- paste(param.txt, " ", sep = "")
  write(param.txt, file = param.file)
  
  # write indfile
  id.order <- sr[[1]]$q.mat$id
  q.mat.df <- do.call(rbind, lapply(sr, function(x) {
    q.mat <- x$q.mat
    rownames(q.mat) <- q.mat$id
    q.mat <- q.mat[id.order, ]
    q.mat <- cbind(row = 1:nrow(q.mat), q.mat)
    rownames(q.mat) <- NULL
    cbind(q.mat[, 1:4], sep = rep(":", nrow(q.mat)), q.mat[, 5:ncol(q.mat)])
  }))
  q.mat.df$pct.miss <- paste("(", q.mat.df$pct.miss, ")", sep = "")
  pop.fac <- factor(q.mat.df$orig.pop)
  pops <- levels(pop.fac)
  q.mat.df$orig.pop <- as.numeric(pop.fac)
  id.fac <- factor(q.mat.df$id)
  ids <- levels(id.fac)
  q.mat.df$id <- as.numeric(id.fac)
  write.table(q.mat.df, file = ind.file, quote = FALSE, sep = " ", 
              row.names = FALSE, col.names = FALSE)
  
  tryCatch({
    err.code <- system(paste("CLUMPP", param.file))
    if(err.code == 127) {
      stop("You do not have CLUMPP installed.")
    } else if(!err.code == 0) {
      stop(paste("Error running CLUMPP. Error code", err.code, "returned."))
    }
  })
  
  if(!file.exists(out.file)) {
    stop(paste("CLUMPP exited without creating output file, '", 
               out.file, "'.", sep = ""))
  }
  out.txt <- scan(out.file, "character", sep = "\n", quiet = TRUE)
  q.mat <- .structureParseQmat(out.txt, pops)
  q.mat$id <- ids[as.numeric(q.mat$id)]
  q.mat$row <- NULL
  
  if(delete.files) {
    files <- c(ind.file, param.file, permutation.file, out.file, misc.file)
    for(f in files) if(file.exists(f)) file.remove(f)
  }
  
  q.mat
}
