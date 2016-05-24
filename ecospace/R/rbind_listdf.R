#' Combine Large Lists of Data Frames.
#'
#' Quickly combine large lists of dataframes into a single data frame by
#' combining them first into smaller data frames. This is only time- and
#' memory-efficient when dealing with large (>100) lists of data frames.
#'
#' @param lists List in which each component is a data frame. Each data frame
#'   must have the same number and types of columns, although the data frames
#'   can have different numbers of rows.
#' @param seq.by Length of small data frames to work with at a time. Default is
#'   100, which timing tests confirm is generally most efficient.
#'
#' @details Rather than combine all list data frames into a single data frame,
#'   this function builds smaller subsets of data frames, and then combines them
#'   together at the end. This process is more time- and memory-efficient.
#'   Timing tests confirm that seq.by=100 is the optimal break-point size. See
#'   examples for confirmation. Function can break large lists into up to
#'   456,976 data frame subparts, giving a warning if requires more subparts.
#'   Only time- and memory-efficient when dealing with large (>100) lists of
#'   data frames.
#'
#' @return Single data frame with number of rows equal to the sum of all data
#'   frames in the list, and columns the same as those of individual list data
#'   frames.
#'
#' @note When called within \code{lapply}, the simulation functions
#'   \code{neutral}, \code{redundancy}, \code{partitioning}, \code{expansion},
#'   and \code{calc_metrics}, which calculates common ecological disparity
#'   (functional diversity) statistics on simulated biotas, produce lists of
#'   data frames. This function is useful for combining these separate lists
#'   into a single data frame for subsequent analyses. This is especially useful
#'   when using the functions within a high-performance computing environment
#'   when submitted as 'embarrassingly parallel' implementations.
#'
#' @seealso \code{\link{neutral}}, \code{\link{redundancy}},
#'   \code{\link{partitioning}}, \code{\link{expansion}}, and
#'   \code{\link{calc_metrics}}
#'
#' @examples
#' nl <- 1000     # List length
#' lists <- vector("list", length=nl)
#' for(i in 1:nl) lists[[i]] <- list(x = rnorm(100), y = rnorm(100))
#' str(lists)
#' object.size(lists)     # ~ 2 MB
#' all <- rbind_listdf(lists)
#' str(all)
#' object.size(all)       # ~ 7 MB
#'
#' # Note each of following can take a few seconds to run
#' # Compare timings:
#' t0 <- Sys.time()
#' all <- rbind_listdf(lists)
#' (Sys.time() - t0)
#'
#' t0 <- Sys.time()
#' all <- rbind_listdf(lists, seq.by=50)
#' (Sys.time() - t0)
#'
#' t0 <- Sys.time()
#' all <- rbind_listdf(lists, seq.by=500)
#' (Sys.time() - t0)
#'
#' # Compare to non-function version
#' all2 <- data.frame()
#' t0 <- Sys.time()
#' for(i in 1:nl) all2 <- rbind(all2, lists[[i]])
#' (Sys.time() - t0)
#'
#' # Build blank ecospace framework to use in simulations
#' ecospace <- create_ecospace(nchar=15, char.state=rep(3, 15), char.type=rep("numeric", 15))
#'
#' # Build 25 samples for neutral model:
#' nreps <- 1:30
#' n.samples <- lapply(X=nreps, FUN=neutral, Sseed=3, Smax=10, ecospace)
#'
#' # Calculate functional diversity metrics for simulated samples
#' n.metrics <- lapply(X=nreps, FUN=calc_metrics, samples=n.samples, Model="neutral", Param="NA")
#' alarm()
#' str(n.metrics)
#'
#' # rbind lists together into a single dataframe
#' all <- rbind_listdf(n.metrics)
#'
#' # Calculate mean dynamics (not reliable given only 10 replicates)
#' means <- n.metrics[[1]]
#' for(n in 1:10) means[n,4:11] <- apply(all[which(all$S==means$S[n]),4:11], 2, mean, na.rm=TRUE)
#' means
#'
#' # Plot statistics as function of species richness, overlaying mean dynamics
#' op <- par()
#' par(mfrow=c(2,4), mar=c(4, 4, 1, .3))
#' attach(all)
#'
#' plot(S, H, type="p", cex=.75, col="gray")
#' lines(means$S, means$H, type="l", lwd=2)
#' plot(S, D, type="p", cex=.75, col="gray")
#' lines(means$S, means$D, type="l", lwd=2)
#' plot(S, M, type="p", cex=.75, col="gray")
#' lines(means$S, means$M, type="l", lwd=2)
#' plot(S, V, type="p", cex=.75, col="gray")
#' lines(means$S, means$V, type="l", lwd=2)
#' plot(S, FRic, type="p", cex=.75, col="gray")
#' lines(means$S, means$FRic, type="l", lwd=2)
#' plot(S, FEve, type="p", cex=.75, col="gray")
#' lines(means$S, means$FEve, type="l", lwd=2)
#' plot(S, FDiv, type="p", cex=.75, col="gray")
#' lines(means$S, means$FDiv, type="l", lwd=2)
#' plot(S, FDis, type="p", cex=.75, col="gray")
#' lines(means$S, means$FDis, type="l", lwd=2)
#'
#' par(op)
#'
#' @export
rbind_listdf <- function(lists=NULL, seq.by=100) {
  nr <- length(lists)
  seq.start <- seq.int(1, nr, by=seq.by)
  lseq <- length(seq.start)
  seq.end <- sort(unique(c((seq.start - 1), nr)))
  seq.end <- seq.end[seq.end >= min(seq.start) & seq.end <= nr]
  seq.start <- seq.start[1:length(seq.end)]
  alphas <- expand.grid(LETTERS[1:26], LETTERS[1:26] , LETTERS[1:26], LETTERS[1:26])
  alphas <- paste(alphas[,4], alphas[,3], alphas[,2], alphas[,1], sep="")
  if(lseq > length(alphas)) stop("only 456,976 temporary variables to store more than that many parts. Make seq.by larger (or modify original rbind_listdf function\n")
  dfs <- paste("df", alphas[seq(lseq)], sep="")
  # rbind as small pieces:
  for (b in 1:lseq) {
    assign(dfs[b], data.frame())# Create blank dfs
    for(c in seq.start[b]:seq.end[b]) {
      assign(dfs[b], rbind(get(dfs[b]), lists[[c]]))
    }
  }
  # rbind all back together:
  out <- data.frame()
  for(b in 1:lseq) { out <- rbind(out, get(dfs[b])) }
  return(out)
}
