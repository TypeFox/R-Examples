#' @title Low Frequency Substitutions
#' @description Check nucleotide sites for low frequency substitutions.
#' 
#' @param x a \code{\link[ape]{DNAbin}} object.
#' @param min.freq minimum frequency of base to be flagged.
#' @param motif.length length of motif around low frequency base to output.
#' @param ... arguments passed from other functions (ignored).
#' 
#' @return data.frame listing id, site number, and motif around low frequency 
#'   base call.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.haps)
#' library(ape)
#' 
#' lowFreqSubs(as.DNAbin(dolph.haps))
#' 
#' @export
#' 
lowFreqSubs <- function(x, min.freq = 3, motif.length = 10, ...) {  
  
  if(!inherits(x, "DNAbin")) stop("'x' must be a DNAbin object")
  x <- as.character(as.matrix(x))
  
  motif.half <- round(motif.length / 2, 0)
  var.sites <- variableSites(x)
  has.min.freq <- apply(var.sites$site.freq, 2, function(site.freq) {
    site.freq <- site.freq[site.freq > 0 & site.freq < min.freq]
    length(site.freq) > 0
  })
  sites.w.min.freq <- var.sites$site.freq[, has.min.freq]
  sites.to.check <- lapply(colnames(sites.w.min.freq), function(site) {
    position <- as.numeric(site)
    site.freq <- sites.w.min.freq[, site]
    bases <- names(site.freq)[which(site.freq > 0 & site.freq < min.freq)]
    site <- x[, position]
    id <- names(site)[which(site %in% bases)]
    to.check <- data.frame(id = id, site = rep(position, length(id)))
    to.check$base <- x[id, position]      
    to.check$motif <- sapply(id, function(i) {
      start.bp <- max(1, position - motif.half)
      end.bp <- min(ncol(x), position + motif.half)
      paste(x[i, start.bp:end.bp], collapse = "")
    })
    to.check
  })
  sites.to.check <- do.call(rbind, sites.to.check)  
  sites.to.check <- sites.to.check[order(sites.to.check$id), ]
  rownames(sites.to.check) <- NULL
  sites.to.check
}