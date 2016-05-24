#' Read in the output of the genosToABH plugin.
#'
#' @param pathToABH The path and filename of the input file.
#' @param nameA Name of the parent represented by "A" in the input file.
#' @param nameB Name of the parent represented by "B" in the input file.
#' @param readPos Should the function attempt read the physical position of
#' markers from the input ?
#'
#' @return A genotype list object which holds the information from the input file.
#'  This list is the fundamental datastructure used by the other functions in this
#'  package. See the vignette for what each item in the list is.
#'
#' @details The input files should be a .csv file holding genotypes as specified by
#'  the qtl package and its "csvs" format.
#'  All characters in the genotype matrix which are not either A,B or H
#'  will be set to N.
#'  If readPos = TRUE (default) marker names must conform to S1_123456 meaning 123456 bp
#'  on chromosome 1. If FALSE, pos is set to NULL and needs to be manually constructed
#'  as shown in the examples. Note that this might throw off some plotting function.
#'
#' @examples \dontrun{genotypes <- readABHgenotypes("./genotypes.csv", "NB", "OL")}
#'
#'  \dontrun{otherGenotypes <- readABHgenotypes("./otherGenotypes.csv", readPos = FALSE)}
#'  #arbitrary position to keep marker order intact
#'  \dontrun{therGenotypes$pos <- 1:length(otherGenotypes$marker_names)}

#' @export
readABHgenotypes <- function(pathToABH, nameA = "A", nameB = "B", readPos = TRUE){

  HapMap <- read.delim(file = pathToABH,
                       stringsAsFactors=FALSE,
                       sep = ",")

  if(readPos == TRUE) {
    pos_temp <- as.integer(sub(pattern = ".+_",
                               x = as.character(colnames(HapMap[,-1])),
                               replacement = "", perl = TRUE)) #replace everything before and including the _ with ""}
  } else {pos_temp <- NULL}

  #get names and factors for each site
  genotypes <- list(
    "ABHmatrix" = as.matrix(HapMap[-1,-1]),
    "chrom" = as.integer(HapMap[1,-1]),
    "individual_names" = as.character(HapMap[-1,1]),
    "marker_names" = as.character(colnames(HapMap[,-1])),
    "pos" = pos_temp,
    "nameA" = nameA,
    "nameB" = nameB
  )

  dimnames(genotypes$ABHmatrix) <- list("individual_names" = HapMap[-1,1],
                                        "marker_names" = colnames(HapMap[,-1]))

  genotypes$ABHmatrix[!(genotypes$ABHmatrix %in% c("A", "B", "H"))] <- "N"

  genotypes
}
