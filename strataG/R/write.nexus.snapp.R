#' @title Write NEXUS File for SNAPP
#' @description Write NEXUS File for SNAPP
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param file the filename the NEXUS file to output.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
write.nexus.snapp <- function(g, file = "snapp.data.nex") {
  if(ploidy(g) != 2) stop("'g' must have diploid data")
  
  num.alleles <- numAlleles(g)
  biallelic <- names(num.alleles)[num.alleles <= 2]
  if(length(biallelic) == 0) {
    warning("No loci are biallelic. No file written.")
    return(NULL)
  }
  
  g <- g[, biallelic, ]
  result <- sapply(locNames(g), function(locus) {
    genotypes <- sapply(indNames(g), function(ind) {
      x <- unlist(as.numeric(loci(g, ids = ind, loci = locus)[, locus]))
      switch(paste(sort(x), collapse = "."), '1.1' = 0, '1.2' = 1, '2.2' = 2, NA)
    })
    as.character(as.numeric(factor(genotypes)) - 1)
  })
  
  result[is.na(result)] <- "?"
  result <- lapply(1:nrow(result), function(i) result[i, ])
  strata <- gsub("[ _]", ".", strata(g))
  id <- gsub("[ _]", ".", indNames(g))
  names(result) <- paste(strata, id, sep = "_")
  
  write.nexus.data(result, file = file)
  
  snapp.file <- scan(file, what = "character", sep = "\n", quiet = TRUE)
  bgn <- grep("BEGIN", snapp.file)
  snapp.file[bgn] <- "BEGIN CHARACTERS;"
  fmt <- grep("FORMAT", snapp.file)
  snapp.file[fmt] <- "  FORMAT DATATYPE=STANDARD MISSING=? GAP=- SYMBOLS=\"012\" LABELS=LEFT TRANSPOSE=NO INTERLEAVE=NO;"
  write(snapp.file, file = file)
  
  invisible(do.call(rbind, result))
}