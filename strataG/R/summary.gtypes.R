#' @name summary,gtypes-method
#' @title Summarize gtypes Object
#' @description Generate a summary of a \code{gtypes} object.
#'  
#' @param object a \linkS4class{gtypes} object.
#' @param x list from summary.gtypes
#' @param ... other arguments (ignored).
#' 
#' @return a list with the following elements:
#' \tabular{ll}{
#'   \code{num.ind} \tab number of individuals.\cr
#'   \code{num.loc} \tab number of loci.\cr
#'   \code{num.strata} \tab number of strata.\cr
#'   \code{allele.freqs} \tab a list with tables of allele frequencies 
#'     by strata.\cr
#'   \code{strata.smry} \tab a by-strata data.frame summarizing haplotypes 
#'     or loci.\cr
#'   \code{locus.smry} \tab a data.frame summarizing each locus for 
#'     non-haploid objects, \code{NULL} for haploid objects.\cr
#'   \code{seq.smry} \tab a summary of the sequence length and base 
#'     frequencies.\cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
setMethod("summary", "gtypes",
  function(object, ...) { 
    x <- object
  
    smry <- list(num.ind = nInd(x), num.loc = nLoc(x), num.strata = nStrata(x))
    smry$allele.freqs <- alleleFreqs(x, by.strata = TRUE)
    
    smry$strata.smry <- t(sapply(strataSplit(x), function(g) {
      c(num.samples = nInd(g),
        num.missing = mean(numMissing(g)),
        num.alleles = mean(numAlleles(g)),
        prop.unique.alleles = mean(propUniqueAlleles(g)),
        heterozygosity = if(ploidy(g) == 1) {
          mean(exptdHet(g))
        } else {
          mean(obsvdHet(g))
        }
      )
    }))
    
    smry$unstratified <- sum(is.na(strata(x))) 
  
    smry$locus.smry <- if(ploidy(x) > 1) summarizeLoci(x) else NULL
    
    smry$seq.smry <- if(!is.null(sequences(x))) {
      sequences <- getSequences(sequences(x), simplify = FALSE)
      do.call(rbind, sapply(sequences, function(dna) {
        dna <- as.matrix(dna)
        dna.len <- unlist(lapply(dna, length))
        len.range <- range(dna.len)
        result <- data.frame(num.seqs = nrow(dna),
                             min.length = len.range[1], 
                             mean.length = round(mean(dna.len)), 
                             max.length = len.range[2]
        )
        cbind(result, rbind(base.freq(dna)))
      }, simplify = FALSE))
    } else NULL
  
    attr(smry, "description") <- x@description
    class(smry) <- c("gtypeSummary", "list")
    smry
})


#' @rdname summary-gtypes-method
#' @export
#' 
print.gtypeSummary <- function(x, ... ) { 
  ind.txt <- paste(x$num.ind, " sample", 
                   ifelse(x$num.ind > 1, "s", ""), sep = "")
  loc.txt <- paste(x$num.loc, " loc", 
                   ifelse(x$num.loc > 1, "i", "us"), sep = "")
  strata.txt <- paste(x$num.strata, " strat", 
                      ifelse(x$num.strata > 1, "a", "um"), sep = "")
  
  cat("\n")
  cat("<<<", attr(x, "description"), ">>>\n")
  cat("\nContents: ")
  cat(ind.txt, loc.txt, strata.txt, sep = ", ")
  cat("\n\nStrata summary:\n")
  print(x$strata.smry)
  if(x$unstratified > 0) cat(x$unstratified, "samples are unstratified\n")
  if(!is.null(x$locus.smry)) {
    cols <- c(1, 3, 5, 7)
    num.rows <- nrow(x$locus.smry)
    if(num.rows > 30) {
      cat("\nLocus summary (first and last 10):\n")
      print(x$locus.smry[1:10, cols, drop = FALSE])
      cat("...\n")
      print(x$locus.smry[(num.rows - 10):num.rows, cols, drop = FALSE])
    } else {
      cat("\nLocus summary:\n")
      print(x$locus.smry[, cols, drop = FALSE])
    }
  }
  if(!is.null(x$seq.smry)) {
    cat("\nSequence summary:\n")
    print(x$seq.smry)
  }
  cat("\n")
  
  invisible(x)
}