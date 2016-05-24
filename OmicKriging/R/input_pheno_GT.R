#' Loads sample phenotype and covariate data into data frame.
#'
#' This function loads a file into a data frame. This file should contain one
#' row per sample in your study, and one column for each covariate and
#' phenotype of interest. Additionally, it requires a header with "IID" for
#' the column of sample IDs, and a unique name for each phenotype and covariate.
#'
#' @param phenoFile File path to the phenotype/covariate file.
#' @param main.pheno Column name of the main phenotype of interest.
#'
#' @return A data frame with dimensions (# of samples) x (# of phenotypes/covar)
#'
#' @keywords input
#'
#' @importFrom utils read.table
#' @export
load_sample_data <- function(phenoFile, main.pheno) {

  ## TODO:: make columns numeric when they can reasonably be converted into numeric

  # load phenotype data
  pheno <- read.table(phenoFile, header=T)
  pheno[main.pheno] <- as.numeric(unlist(pheno[main.pheno]))
  rownames(pheno) <- pheno$IID

  return(pheno)
}

