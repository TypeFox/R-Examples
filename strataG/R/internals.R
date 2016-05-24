#' @importFrom parallel detectCores makePSOCKcluster makeForkCluster
#' 
.setupClusters <- function(num.cores = NULL) {
  # setup clusters
  if(is.null(num.cores)) num.cores <- detectCores() - 1
  if(is.na(num.cores)) num.cores <- 1
  num.cores <- max(1, num.cores)
  num.cores <- min(num.cores, detectCores() - 1)
  if(num.cores > 1) {
    is.windows <- .Platform$OS.type == "windows"
    cl.func <- ifelse(is.windows, makePSOCKcluster, makeForkCluster)
    cl.func(num.cores)
  } else NULL
}


#' @importFrom utils combn
#' 
.strataPairs <- function(g) {
  if(nlevels(strata(g)) < 2) return(NULL)
  strata.vec <- sort(unique(as.character(strata(g))))
  strata.pairs <- t(combn(strata.vec, 2))
  colnames(strata.pairs) <- c("strata.1", "strata.2")
  as.data.frame(strata.pairs, stringsAsFactors = FALSE)
}

