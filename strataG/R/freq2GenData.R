#' @title Convert Haplotype Frequency Matrices
#' @description Create a data.frame of stratified individuals and their 
#'   haplotypes from a frequency table
#' 
#' @param freq.mat a matrix or data.frame containing haplotypic 
#'   frequencies with strata as column names.
#' @param hap.col a number giving the column providing haplotype 
#'   labels or a vector the same length as freq.mat.
#' @param freq.col a number giving the first column containing 
#'   haplotype frequencies.
#' @param id.label character to label sample IDs with in resulting 
#'   data.frame.
#' @param hap.label character to label haplotypes with in 
#'   resulting data.frame.
#' 
#' @return a data.frame with one row per sample and columns for id, strata, 
#'   and haplotype.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' hap.freqs <- data.frame(
#'   haps = c("hap1", "hap2", "hap3"),
#'   pop1 = rmultinom(1, 50, prob = c(0.1, 0.2, 0.7)),
#'   pop2 = rmultinom(1, 25, prob = c(0.5, 0.4, 0.1))
#' )
#' 
#' gen.data <- freq2GenData(hap.freqs, hap.col = 1, freq.col = 2)
#' 
#' x <- df2gtypes(gen.data, ploidy = 1)
#' summary(x)
#' 
#' @importFrom stats complete.cases
#' @export
#' 
freq2GenData <- function(freq.mat, hap.col, freq.col, id.label = NULL, 
                         hap.label = NULL) {
  haps <- if(length(hap.col) == 1) freq.mat[, hap.col] else hap.col
  df <- do.call(rbind, lapply(freq.col:ncol(freq.mat), function(i) {
    freq.sub <- cbind(haps, freq.mat[, i])
    freq.sub <- freq.sub[complete.cases(freq.sub), , drop = FALSE]
    haps <- rep(freq.sub[, 1], times = as.numeric(freq.sub[, 2]))
    if(!is.null(hap.label)) haps <- paste(hap.label, haps, sep = ".")   
    data.frame(strata = colnames(freq.mat)[i], haplotype = haps, 
               stringsAsFactors = FALSE)
  }))
  ids <- 1:nrow(df)
  if(!is.null(id.label)) ids <- paste(id.label, ids, sep = ".")
  rownames(df) <- ids
  cbind(id = ids, df, stringsAsFactors = FALSE)
}