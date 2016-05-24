#' Calculate multilocus heterozygosity (MLH) 
#' 
#' sMLH is defined as the total number of heterozygous loci in an individual divided 
#' by the sum of average observed heterozygosities in the population over the 
#' subset of loci successfully typed in the focal individual.
#'
#' @param genotypes data.frame with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote), 1 (heterozygote) and \code{NA} (missing)
#'
#' @return
#' Vector of individual standardized multilocus heterozygosities
#'
#' @references
#' Coltman, D. W., Pilkington, J. G., Smith, J. A., & Pemberton, J. M. (1999). 
#' Parasite-mediated selection against inbred Soay sheep in a free-living, 
#' island population. Evolution, 1259-1267.
#'
#' @author Martin A. Stoffel (martin.adam.stoffel@@gmail.com) 
#'        
#' @examples
#' data(mouse_msats)
#' genotypes <- convert_raw(mouse_msats)
#' het <- sMLH(genotypes)
#'
#' @export
#'

sMLH <- function(genotypes) {
    # transform to matrix
    genes <- as.matrix(genotypes)
    # genes[is.na(genes)] <- -1
    # get logical matrix of non-missing values as TRUE
    typed <- (genes == 1) | (genes == 0)
    typed[is.na(typed)] <- FALSE
    # initialise
    nloc <- ncol(genes)
    nind <- nrow(genes)
    typed_sum <- colSums(typed, na.rm = TRUE)
    # heterozygosity per locus
    het_loc <- colSums(genes == 1, na.rm = TRUE) / typed_sum
    # replicate vector to matrix
    het_loc_mat <- matrix(het_loc, nrow = nind, ncol = nloc, byrow = TRUE)
    het_loc_mat[!typed] <- 0
    mh <- rowSums(het_loc_mat, na.rm = TRUE)
    N  <- rowSums(typed, na.rm = TRUE)
    H  <- rowSums(genes == 1, na.rm = TRUE)
    sMLH <- (H/N)/(mh/N)
    names(sMLH) <- row.names(genotypes)
    sMLH
}


#' Calculate multilocus heterozygosity (MLH) 
#' 
#' MLH is defined as the total number of heterozygous loci in an individual divided 
#' by the number of loci typed in the focal individual. An MLH of 0.5 thus means
#' that 50 percent of an indiviudals loci are heterozygous.
#'
#' @param genotypes data.frame with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote), 1 (heterozygote) and NA (missing)
#'
#' @return
#' Vector of individual multilocus heterozygosities
#'
#' @references
#' Coltman, D. W., Pilkington, J. G., Smith, J. A., & Pemberton, J. M. (1999). 
#' Parasite-mediated selection against inbred Soay sheep in a free-living, 
#' island population. Evolution, 1259-1267.
#'
#' @author Martin A. Stoffel (martin.adam.stoffel@@gmail.com) 
#' 
#' @examples
#' data(mouse_msats)
#' genotypes <- convert_raw(mouse_msats)
#' het <- MLH(genotypes)
#' 
#' @export
#' 
MLH <- function(genotypes) {
    # transform to matrix
    genes <- as.matrix(genotypes)
    # genes[is.na(genes)] <- -1
    # get logical matrix of non-missing values as TRUE
    typed <- (genes == 1) | (genes == 0)
    typed[is.na(typed)] <- FALSE
    # initialise
    nloc <- ncol(genes)
    nind <- nrow(genes)
    typed_sum <- colSums(typed, na.rm = TRUE)
    # heterozygosity per locus
    het_loc <- colSums(genes == 1, na.rm = TRUE) / typed_sum
    # replicate vector to matrix
    het_loc_mat <- matrix(het_loc, nrow = nind, ncol = nloc, byrow = TRUE)
    het_loc_mat[!typed] <- 0
    mh <- rowSums(het_loc_mat, na.rm = TRUE)
    N  <- rowSums(typed, na.rm = TRUE)
    H  <- rowSums(genes == 1, na.rm = TRUE)
    MLH <- (H/N)
}

