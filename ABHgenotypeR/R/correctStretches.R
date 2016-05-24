#' Correct short miscalled stretches based on flanking alleles.
#'
#' @param inputGenos A genotypes list object.
#' @param maxHapLength The maximum length of stretches flanked
#'   by non-heterzygous sites that are changed. If set to 1
#'   (default) only AXA or BXB will be corrected. If set to 2, both AXA and AXYA
#'   (or BXB and BXYB) will be corrected.
#'
#' @return A genotype object in which short miscalled stretches are
#'   corrected if both flanking alleles match.
#'
#' @examples \dontrun{corrStretchGenos <- correctStretches(genotypes, maxHapLength = 3)}

#' @export
correctStretches <- function(inputGenos = "genotypes",
                             maxHapLength = 1) {

  geno_raw <- inputGenos$ABHmatrix
  #setup a matrix for correcting errors
  geno_correctedErr <- matrix(0,
                              nrow = nrow(geno_raw),
                              ncol = ncol(geno_raw))
  #make reg expressions for A
  patExprA <- NULL # a character vector which holds regexp
  for(i in 1:maxHapLength) {
    patExprA[i] <- paste("(A)([BHN]{",i,"})(?=A)", sep = "")
  }
  replExprA <- NULL # a character vector which holds regexp
  for(i in 1:maxHapLength) {
    replExprA[i] <- paste("\\1",
                          paste(rep("\\1",i), sep = "", collapse = ""),
                          sep = "", collapse = "")
  }

  #make reg expressions for B
  patExprB <- NULL # a character vector which holds regexp
  for(i in 1:maxHapLength) {
    patExprB[i] <- paste("(B)([AHN]{",i,"})(?=B)", sep = "")
  }
  replExprB <- NULL # a character vector which holds regexp
  for(i in 1:maxHapLength) {
    replExprB[i] <- paste("\\1",
                          paste(rep("\\1",i), sep = "", collapse = ""),
                          sep = "", collapse = "")
  }

  #in geno_raw replace errors
  for(chrom_count in unique(inputGenos$chrom)) {
    geno_temp <- geno_raw[,inputGenos$chrom == chrom_count]
    for (row_count in 1:nrow(geno_correctedErr)) {
      for(HapLen in 1:length(patExprA)) { #replace with A
        if(HapLen == 1) { #first replace in genotemp, then in geno_correc
          geno_correctedErr[row_count, inputGenos$chrom == chrom_count] <- substring(gsub(x = paste(geno_temp[row_count,],
                                                                                                    collapse = ""),
                                                                                          pattern = patExprA[HapLen],
                                                                                          replacement = replExprA[HapLen],
                                                                                          perl = TRUE),
                                                                                     1:ncol(geno_temp),
                                                                                     1:ncol(geno_temp))
        } else {
          geno_correctedErr[row_count, inputGenos$chrom == chrom_count] <- substring(gsub(x = paste(geno_correctedErr[row_count, inputGenos$chrom == chrom_count],
                                                                                                    collapse = ""),
                                                                                          pattern = patExprA[HapLen],
                                                                                          replacement = replExprA[HapLen],
                                                                                          perl = TRUE),
                                                                                     1:ncol(geno_temp),
                                                                                     1:ncol(geno_temp))
        }
      }
      for(HapLen in 1:length(patExprB)) { #replace with B, only in geno_correctedErr
        geno_correctedErr[row_count, inputGenos$chrom == chrom_count] <- substring(gsub(x = paste(geno_correctedErr[row_count, inputGenos$chrom == chrom_count],
                                                                                                  collapse = ""),
                                                                                        pattern = patExprB[HapLen],
                                                                                        replacement = replExprB[HapLen],
                                                                                        perl = TRUE),
                                                                                   1:ncol(geno_temp),
                                                                                   1:ncol(geno_temp))
      }
    }
  }

  outputGenos <-inputGenos
  outputGenos$ABHmatrix <- geno_correctedErr
  dimnames(outputGenos$ABHmatrix) <- list("individual_names" = inputGenos$individual_names,
                                          "marker_names" = inputGenos$marker_names)

  reportGenos(inputGenos)
  cat(paste("\n"))
  reportGenos(outputGenos)

  outputGenos
}

