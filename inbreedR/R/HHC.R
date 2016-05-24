#' Calculates heterzygosity-heterozygosity correlations with 
#' standardized multilocus heterozygosities (sMLH)
#' 
#' Loci are randomly devided into two equal groups and the correlation coefficient
#' between the resulting sMLH values is calculated.
#'
#' @param genotypes \code{data.frame} with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote), 1 (heterozygote) and \code{NA} (missing)
#' @param reps number of repetitions, i.e. splittings of the dataset 
#' @param CI size of the confidence interval around the mean het-het correlation (default is 0.95)
#'
#' @return 
#' \item{call}{function call.}
#' \item{HHC_vals}{vector of HHC's obtained by randomly splitting the dataset}
#' \item{summary_exp_r2}{r2 mean and sd for each number of subsetted loci}
#' \item{nobs}{number of observations}
#' \item{nloc}{number of markers}
#' 
#' @references
#' Balloux, F., Amos, W., & Coulson, T. (2004). Does heterozygosity estimate inbreeding
#' in real populations?. Molecular Ecology, 13(10), 3021-3031.
#' 
#' @author Martin A. Stoffel (martin.adam.stoffel@@gmail.com) 
#'        
#' @examples
#' data(mouse_msats)
#' genotypes <- convert_raw(mouse_msats)
#' (out <- HHC(genotypes, reps = 100, CI = 0.95))
#'
#' @export
#'
#'

HHC <- function(genotypes, reps = 100, CI = 0.95) {
    # transform to matrix
    genotypes <- as.matrix(genotypes)
    # initialise
    loci <- 1:ncol(genotypes)
    loc_num <- ncol(genotypes)
    
    calc_cor <- function(num_iter, genotypes) {
        new_ord <- sample(loci)
        sMLH1 <- sMLH(genotypes[, new_ord[1:floor(loc_num/2)]])
        sMLH2 <- sMLH(genotypes[, new_ord[(floor(loc_num/2) + 1):loc_num]])
        het_het_cor <- stats::cor(sMLH1, sMLH2)
        
        if (num_iter == 1) {
            cat("\n", "starting het-het correlations")
        } else if (num_iter == reps) {
            cat("\n", "done")
        } else if (num_iter %% 5 == 0) {
            cat("\n", num_iter, "iterations done")
        }
        return(het_het_cor)
    }
    HHC_vals <- vapply(1:reps, calc_cor, numeric(1), genotypes)
    
    CI_HHC <- stats::quantile(HHC_vals, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)

    res <- list(call = match.call(),
                HHC_vals = HHC_vals,
                CI_HHC = CI_HHC,
                nobs = nrow(genotypes), 
                nloc = ncol(genotypes))
    
    class(res) <- "inbreed"
    return(res)
}




