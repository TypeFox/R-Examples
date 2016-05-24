#' Impute missing genotypes based on flanking alleles
#'
#' @param inputGenos A genotypes list object.
#'
#' @return A genotype object in which missing data is imputed based on flanking
#'   alleles. Any number of N is replaced by either A, B or N if the alleles which flank the N match
#'
#' @examples \dontrun{imputedGenos <- imputeByFlanks(genotypes)}

#' @export
imputeByFlanks <- function (inputGenos = "genotypes") {

  geno_raw <- inputGenos$ABHmatrix
  #setup empty matrix for imputed genotypes
  geno_imp <- matrix(0,
                     nrow = nrow(geno_raw),
                     ncol = ncol(geno_raw))

  #loop through chromosome
  for(chrom_count in unique(inputGenos$chrom)) {
    geno_temp <- geno_raw[,inputGenos$chrom == chrom_count]
    #loop through rows, replace Ns flanked by the same parent/het
    for (row_count in 1:nrow(geno_temp)){
      geno_imp[row_count, inputGenos$chrom == chrom_count] <- substring(gsub(x = paste(geno_temp[row_count,],
                                                                            collapse = ""),
                                                                  pattern = "(?:A|(?<!^)\\G)\\KN(?=N*A)",
                                                                  replacement = "A",
                                                                  perl = TRUE),
                                                             1:ncol(geno_temp),
                                                             1:ncol(geno_temp))
      geno_imp[row_count, inputGenos$chrom == chrom_count] <- substring(gsub(x = paste(geno_imp[row_count, inputGenos$chrom == chrom_count],
                                                                            collapse = ""),
                                                                  pattern = "(?:B|(?<!^)\\G)\\KN(?=N*B)",
                                                                  replacement = "B",
                                                                  perl = TRUE),
                                                             1:ncol(geno_temp),
                                                             1:ncol(geno_temp))
      geno_imp[row_count, inputGenos$chrom == chrom_count] <- substring(gsub(x = paste(geno_imp[row_count, inputGenos$chrom == chrom_count],
                                                                            collapse = ""),
                                                                  pattern = "(?:H|(?<!^)\\G)\\KN(?=N*H)",
                                                                  replacement = "H",
                                                                  perl = TRUE),
                                                             1:ncol(geno_temp),
                                                             1:ncol(geno_temp))
    }
  }
  outputGenos <-inputGenos
  outputGenos$ABHmatrix <- geno_imp
  dimnames(outputGenos$ABHmatrix) <- list("individual_names" = inputGenos$individual_names,
                                        "marker_names" = inputGenos$marker_names)
  reportGenos(inputGenos)
  cat(paste("\n"))
  reportGenos(outputGenos)

  outputGenos
}
