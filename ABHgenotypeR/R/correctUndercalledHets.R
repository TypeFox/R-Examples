#' Correct undercalled heterozygous sites based on flanking alleles.
#'
#' @param inputGenos A genotypes list object.
#' @param maxHapLength The maximum length of not heterozygous stretches flanked
#'   by heterzygous sites that are changed to heterozygous. If set to 1
#'   (default) only HAH or HBH will be corrected. If set to 2, both HAH and HAAH
#'   (or HBH and HBBH) will be corrected.
#'
#' @return A genotype object in which undercalled heterozygous sites are
#'   corrected if both flanking alleles match.
#'
#' @examples \dontrun{corrUndHetsGenos <- correctUndercalledHets(genotypes, maxHapLength = 3)}

#' @export
correctUndercalledHets <- function(inputGenos = "genotypes",
                                   maxHapLength = 1) {

  geno_raw <- inputGenos$ABHmatrix
  #setup a matrix for correcting undercalled hets
  geno_correctedHets <- matrix(0,
                               nrow = nrow(geno_raw),
                               ncol = ncol(geno_raw))
  #make reg expressions for H
  patExprH <- NULL # a character vector which holds regexp
  for(i in 1:maxHapLength) {
    patExprH[i] <- paste("(H)([ABN]{",i,"})(?=H)", sep = "")
  }
  replExprH <- NULL # a character vector which holds regexp
  for(i in 1:maxHapLength) {
    replExprH[i] <- paste("\\1",
                          paste(rep("\\1",i), sep = "", collapse = ""),
                          sep = "", collapse = "")
  }


  #in geno_raw replace undercalled hets
  for(chrom_count in unique(inputGenos$chrom)) {
    geno_temp <- geno_raw[,inputGenos$chrom == chrom_count]
    for (row_count in 1:nrow(geno_correctedHets)) { #first replace in genotemp, then in geno_corrected
      for(HapLen in 1:length(patExprH)) {
        if(HapLen == 1) {
          geno_correctedHets[row_count, inputGenos$chrom == chrom_count] <- substring(gsub(x = paste(geno_temp[row_count,],
                                                                                          collapse = ""),
                                                                                pattern = patExprH[HapLen],
                                                                                replacement = replExprH[HapLen],
                                                                                perl = TRUE),
                                                                           1:ncol(geno_temp),
                                                                           1:ncol(geno_temp))
        } else {
          geno_correctedHets[row_count, inputGenos$chrom == chrom_count] <- substring(gsub(x = paste(geno_correctedHets[row_count, inputGenos$chrom == chrom_count],
                                                                                          collapse = ""),
                                                                                pattern = patExprH[HapLen],
                                                                                replacement = replExprH[HapLen],
                                                                                perl = TRUE),
                                                                           1:ncol(geno_temp),
                                                                           1:ncol(geno_temp))
        }
      }
    }
  }
#put feedbackfunction here
  outputGenos <-inputGenos
  outputGenos$ABHmatrix <- geno_correctedHets
  dimnames(outputGenos$ABHmatrix) <- list("individual_names" = inputGenos$individual_names,
                                       "marker_names" = inputGenos$marker_names)

  reportGenos(inputGenos)
  cat(paste("\n"))
  reportGenos(outputGenos)

  outputGenos
  }

