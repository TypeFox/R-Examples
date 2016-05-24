#' @title Quality Assurance/Quality Control
#' @description Conducts a suite of QA/QC tests. Summarizes missing data and 
#'   homozygosity by individual and locus, and looks for duplicate genotypes 
#'   (see \code{\link{dupGenotypes}}). For sequence data, identifies low 
#'   frequency substitutions (see \code{\link{lowFreqSubs}}), and computes 
#'   haplotype likelihoods (see \code{\link{haplotypeLikelihoods}}).
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param label optional label for output folder and prefix for files.
#' @param ... optional arguments to pass on to summary functions.
#' 
#' @return Files are written for by-sample and by-locus summaries, and duplicate 
#'   genotypes if any are found. If sequences are present, files are written 
#'   identifying low frequency substitutions and haplotype likelihoods.\cr
#'   The return value is a list with the following elements:
#'   
#'   \tabular{ll}{
#'      \code{by.sample} \tab data.frame of by-sample summaries.\cr
#'      \code{by.locus} \tab data.frame of by-locus summaries.\cr
#'      \code{dup.df} \tab data.frame identifying potential duplicates.\cr
#'      \code{by.seq} \tab list of low frequency substitutions and haplotype 
#'        likelihoods for each gene.\cr
#'    }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom utils write.csv
#' @export
#' 
qaqc <- function(g, label = NULL, ...) {
  cat("\n")
  cat(format(Sys.time()), ": Individual summaries\n")
  num.loc <- nLoc(g)
  by.sample <- do.call(rbind, lapply(indNames(g), function(id) {
    smry <- sapply(locNames(g), function(loc) {
      gt <- loci(g, id, loc)
      missing <- any(is.na(gt))
      hmzgs <- if(missing) NA else length(unique(unlist(gt))) == 1
      c(missing = missing, hmzgs = hmzgs)
    })
    
    missing <- sum(smry["missing", ], na.rm = TRUE)
    data.frame(sample = id, 
      num.loci.missing.genotypes = missing,
      pct.loci.missing.genotypes = missing / num.loc,
      pct.loci.homozygous = mean(smry["hmzgs", ], na.rm = TRUE)
    )
  }))
  
  cat(format(Sys.time()), ": Locus summaries\n")
  by.locus <- summarizeLoci(g, TRUE)
  by.locus$All <- summarizeLoci(g)
  
  cat(format(Sys.time()), ": Duplicate genotypes\n")
  dup.df <- dupGenotypes(g, ...)
  
  # Sequence summaries
  by.seq <- if(!is.null(sequences(g))) {
    cat(format(Sys.time()), ": Sequence summaries\n")
    sapply(locNames(g), function(x) {
      x.seqs <- getSequences(sequences(g), loci = x, simplify = TRUE)
      list(
        low.freq.subs = lowFreqSubs(x.seqs, ...),
        hap.likelihoods = haplotypeLikelihoods(x.seqs, ...)
      )
    }, simplify = FALSE)
  } else NULL
      
  cat(format(Sys.time()), ": Writing files\n")
  # Write summaries to files
  label <- if(is.null(label)) description(g) else label
  label <- gsub("[[:punct:]]", ".", label)
  fname <- paste(label, ".sample.summary.csv", sep = "")
  write.csv(by.sample, file = fname, row.names = FALSE)

  for(x in names(by.locus)) {
    fname <- paste(label, "locus.summary", x, "csv", sep = ".")
    by.locus[[x]] <- data.frame(
      locus = rownames(by.locus[[x]]), by.locus[[x]], stringsAsFactors = FALSE
    )
    write.csv(by.locus[[x]], file = fname, row.names = FALSE)
  }  
  
  if(!is.null(dup.df)) {
    fname <- paste(label, ".duplicate.samples.csv", sep = "")
    write.csv(dup.df, file = fname, row.names = FALSE)
  }
  
  if(!is.null(by.seq)) {
    for(x in names(by.seq)) {
      fname <- paste(label, "low.freq.subs", x, "csv", sep = ".")
      write.csv(by.seq[[x]]$low.freq.subs, file = fname, row.names = FALSE)
      fname <- paste(label, "haplotype.likelihoods", x, "csv", sep = ".")
      hl <- cbind(delta.logLik = by.seq[[x]]$hap.likelihoods)
      write.csv(hl, file = fname)
    }
  }
  
  cat("\n")
  invisible(list(by.sample = by.sample, by.locus = by.locus, dup.df = dup.df,
                 by.seq = by.seq))
}