#' @title Fixed Differences
#' @description Summarize fixed base pair differences between strata.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param count.indels logical. Count indels when evaluating sites for fixed 
#'   differences?
#' @param consec.indels.as.one logical. If \code{count.indels} is TRUE, count 
#'   consecutive indels as a a single indel?
#' @param bases a character vector of valid bases to consider.
#' 
#' @return a list with components:
#' \tabular{ll}{
#'   \code{sites} \tab list of sites with fixed differences for each pair of 
#'     strata \cr
#'   \code{num.fixed} \tab data.frame of number of sites fixed between each 
#'     pair of strata \cr
#' }
#' 
#' @author Eric Archer <eric.archer@@noaa.gov>
#' 
#' @seealso \code{\link{fixedSites}}, \code{\link{variableSites}}
#' 
#' @export
#' 
fixedDifferences <- function(g, count.indels = TRUE, 
                             consec.indels.as.one = TRUE, 
                             bases = c("a", "c", "g", "t", "-")) {  
  # get fixed sites for each strata
  strata.gtypes <- strataSplit(g, remove.sequences = TRUE)
  fixed.sites <- sapply(strata.gtypes, function(strata) {
    fixedSites(strata, bases)
  }, simplify = FALSE)
  
  # for each pair of strata, return matrix of sites with fixed differences
  # fixed differences are sites which aren't variable in a strata and 
  #   different between strata
  strata.pairs <- as.matrix(.strataPairs(g))
  pair.fixed.diff <- lapply(1:nrow(strata.pairs), function(p) {
    sp <- strata.pairs[p, ]
    fixed.1 <- fixed.sites[[sp[1]]]
    fixed.2 <- fixed.sites[[sp[2]]]
    shared.sites <- intersect(names(fixed.1), names(fixed.2))
    seq.mat <- rbind(fixed.1[shared.sites], fixed.2[shared.sites])
    rownames(seq.mat) <- sp
    if(ncol(seq.mat) == 0) return(list(sites = seq.mat, num.fixed = 0))
    
    # count fixed nucleotides (don't count sites with indels)
    nucs.diff <- apply(seq.mat, 2, function(x) {
      x[1] != x[2] & !any(x %in% c("-", "."))
    })
    sites <- which(nucs.diff)
    num.fixed <- sum(nucs.diff)
    
    # count fixed indels
    if(count.indels) {
      is.diff.indel <- apply(seq.mat, 2, function(x) {
        any(x == "-") & x[1] != x[2]
      })
      indels.fixed <- if(consec.indels.as.one) {
        # mark non-shared indels as '*'
        indel.diff <- apply(seq.mat, 2, function(i) {
          has.indel <- i == "-"
          if(any(has.indel) & i[1] != i[2]) i[has.indel] <- "*"
          i
        })
        # count runs of non-shared indels in either sequence
        f1.rle <- rle(indel.diff[1, ])    
        f2.rle <- rle(indel.diff[2, ])    
        sum(f1.rle$values == "*") + sum(f2.rle$values == "*")
      } else sum(is.diff.indel)
      sites <- sort(c(sites, which(is.diff.indel)))
      num.fixed <- sum(num.fixed, indels.fixed)
    }  
    list(sites = seq.mat[, sites, drop = F], num.fixed = num.fixed)
  })
  
  # compile sites with fixed differences
  sites <- lapply(pair.fixed.diff, function(x) x$sites)
  names(sites) <- apply(strata.pairs, 1, paste, collapse = " v. ")
  
  # count number of fixed differencs between pairs
  num.fixed <- sapply(pair.fixed.diff, function(x) x$num.fixed)
  strata.pairs <- cbind(data.frame(strata.pairs), num.fixed = num.fixed)
  
  list(sites = sites, num.fixed = strata.pairs)
}