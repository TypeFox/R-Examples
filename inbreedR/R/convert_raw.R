#' Genotype format converter
#'
#' Turns raw genotype data into 0 (homozygote), 1 (heterozygote) and \code{NA} (missing), which is the working format for 
#' the \code{inbreedR} functions.
#' A raw genotype matrix has individuals in rows and each locus in two adjacent columns. Individual ID's can be rownames.
#' Type data(mouse_msats) for an example raw genotype data frame.
#'
#' @param genotypes Raw genotype \code{data.frame} or \code{matrix}. Rows represent individuals and each locus has two adjacent columns. 
#'        Alleles within loci can be coded as numbers (e.g. microsatellite length) or characters (e.g. "A", "T")
#'        See data(mouse_msat) for an example. Missing values should be coded as \code{NA}.
#'
#' @return \code{data.frame} object with 0 (homozygote), 1 (heterozygote) and \code{NA} (missing data).
#'         Each locus is a column and each individual is a row.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'
#' @examples
#' # Mouse microsatellite data with missing values coded as NA
#' data(mouse_msats)
#' genotypes <- convert_raw(mouse_msats)
#' head(genotypes)
#'
#' @export


convert_raw <- function(genotypes) {

check_het <- function(x) {

        if ((length(x) %% 2) == 1) {
                s1 <- seq(1, (length(x)-1), 2)
                newx <- as.vector(rep(NA, ceiling(length(x)/2)))
                }
        if ((length(x) %% 2) == 0) {
                s1 <- seq(1, length(x), 2)
                newx <- as.vector(rep(NA, length(x)/2))
                }

        count <- 1
        for(i in s1){
                if ((is.na(x[i]) | is.na(x[i + 1]))) {
                        newx[count] <- NA
                        count <- count + 1
                } else if (x[i] == x[i + 1]) {
                        newx[count] <- 0
                        count <- count + 1
                } else if (x[i] != x[i + 1]) {
                        newx[count] <- 1
                        count <- count + 1
                }
        }
        newx
}

# original full data matrix
origin <- as.data.frame(t(apply(genotypes, 1, check_het)))
row.names(origin) <- row.names(genotypes)
origin
}
