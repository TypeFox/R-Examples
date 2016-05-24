#' @title Nucleotide Divergence
#' @description Calculate Nei's dA between strata, and distributions of 
#'   between- and within-strata nucleotide divergence (sequence distance).
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param probs a numeric vector of probabilities of the pairwise distance 
#'   distributions with values in \code{0:1}.
#' @param model evolutionary model to be used. see \code{\link[ape]{dist.dna}} 
#'   for options.
#' @param ... other arguments passed to \code{\link[ape]{dist.dna}}.
#' 
#' @return a list with summaries of the \code{within} and \code{between} strata 
#'   pairwise distances including Nei's dA. 
#'   
#' @references Nei, M., and S. Kumar (2000) Molecular Evolution and 
#'   Phylogenetics. Oxford University Press, Oxford. (dA: pp. 256, eqn 12.67)
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.seqs)
#' strata <- dolph.strata$fine
#' names(strata) <- dolph.strata$ids
#' dloop <- sequence2gtypes(dolph.seqs, strata, seq.names = "dLoop")
#' 
#' nucleotideDivergence(dloop)
#' 
#' @aliases dA
#' @importFrom stats quantile
#' @export
#' 
nucleotideDivergence <- function(g, probs = c(0, 0.025, 0.5, 0.975, 1), model = "raw", ...) { 
  if(ploidy(g) > 1) stop("'g' must be haploid")
  if(is.null(g@sequences)) stop("'g' must have sequences")
  
  pair.dist.summary <- function(haps, d, probs) {
    pws.dist <- apply(haps, 1, function(h) {
      if(any(is.na(h))) return(NA)
      d[h[1], h[2]]
    })
    dist.quant <- quantile(pws.dist, probs, na.rm = TRUE)
    names(dist.quant) <- paste("pct.", probs, sep = "")
    c(mean = mean(pws.dist, na.rm = TRUE), dist.quant)
  }
  
  g <- g[, , strataNames(g)]
  st.pairs <- .strataPairs(g)
  if(!is.null(st.pairs)) st.pairs <- as.matrix(st.pairs)
  st <- strata(g)
  hap.dist <- lapply(
    getSequences(sequences(g), simplify = FALSE), 
    dist.dna, model = model, as.matrix = TRUE, ...
  )
  
  result <- lapply(locNames(g), function(loc) {
    within.dist <- do.call(rbind, tapply(names(st), st, function(ids) {
      haps <- as.character(loci(g, ids = ids, loci = loc)[, 1])
      pair.dist.summary(t(combn(haps, 2)), hap.dist[[loc]], probs)
    }))
  
    between.dist <- if(is.null(st.pairs)) NA else {
      btwn <- t(apply(st.pairs, 1, function(sp) {
        ids1 <- names(st)[which(st == sp[1])]
        ids2 <- names(st)[which(st == sp[2])]
        h1 <- as.character(loci(g, ids = ids1, loci = loc)[, 1])
        h2 <- as.character(loci(g, ids = ids2, loci = loc)[, 1])
        btwn <- pair.dist.summary(expand.grid(h1, h2), hap.dist[[loc]], probs)
        dA <- btwn["mean"] - (sum(within.dist[sp, "mean"], na.rm = TRUE) / 2)
        c(dA = unname(dA), btwn)
      }))
      data.frame(st.pairs, btwn, stringsAsFactors = FALSE)
    }
    
    list(within = within.dist, between = between.dist) 
  })
  names(result) <- locNames(g)
  result
}