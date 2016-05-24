#' Estimating g2 from microsatellite data
#'
#' @param genotypes \code{data.frame} with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote), 1 (heterozygote) and \code{NA} (missing)
#' @param nperm Number of permutations for testing the hypothesis that the empirical g2-value is higher than the g2 for random associations between 
#'        individuals and genotypes.
#' @param nboot Number of bootstraps for estimating a confidence interval
#' @param CI Confidence interval (default to 0.95)
#' @param verbose If FALSE, nothing will be printed to show the status of bootstraps and permutations.
#'
#' @details Calculates g2 from smaller datasets. The underlying formula is compationally expensive 
#'          due to double summations over all paits of loci (see David et al. 2007). 
#'          Use convert_raw to convert raw genotypes (with 2 columns per locus) into
#'          the required format.
#'          
#' @return
#' g2_microsats returns an object of class "inbreed".
#' The functions `print` and `plot` are used to print a summary and to plot the distribution of bootstrapped g2 values and CI.
#' 
#' An `inbreed` object from \code{g2_microsats} is a list containing the following components:
#' 
#' \item{call}{function call.}
#' \item{g2}{g2 value}
#' \item{p_val}{p value from permutation test}
#' \item{g2_permut}{g2 values from permuted genotypes}
#' \item{g2_boot}{g2 values from bootstrap samples}
#' \item{CI_boot}{confidence interval from bootstraps}
#' \item{se_boot}{standard error of g2 from bootstraps}
#' \item{nobs}{number of observations}
#' \item{nloc}{number of markers}
#'
#' @references
#' David, P., Pujol, B., Viard, F., Castella, V. and Goudet, J. (2007),
#' Reliable selfing rate estimates from imperfect population genetic data. Molecular Ecology,
#' 16: 2474
#'
#' @author Martin A. Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Mareike Esser (messer@@techfak.uni-bielefeld.de)
#'
#'
#' @examples
#' data(mouse_msats)
#' # tranform raw genotypes into 0/1 format
#' genotypes <- convert_raw(mouse_msats)
#' (g2_mouse <- g2_microsats(genotypes, nperm = 1000, nboot = 100, CI = 0.95))
#'
#' @export
#'


g2_microsats <- function(genotypes, nperm = 0, nboot = 0, CI = 0.95, verbose = TRUE) {
    
    genotypes <- as.matrix(genotypes)
    # transpose
    origin <- (t(genotypes)) 
    origin[(origin!=0) & (origin!=1)] <- -1
    origin[is.na(origin)] <- -1
    
    calc_g2 <- function(origin, perm = 1, boot = 1) {
        # define matrix with 1 for missing and 0 for all others
        m <- origin
        m[m == 1] <- 0
        m[m == -1] <- 1
        # H matrix with 0 for -1
        h <- origin
        h[(h!=0)&(h!=1)] <- 0
        
        n <- ncol(origin) # number of individuals
        l <- nrow(origin) # number of loci
        
        # mij: individuals missing on i and j ?s locus
        m_ij <- (m %*% t(m))
        # missings per locus
        m_loc <- rowSums(m)
        
        # numerator --------------------------------------------------------------------
        # pij entry says total amount of individuals that are heterozygous at locus i and locus j
        p <- h %*% t(h)
        numerator_mat <- matrix(rep(0, l*l), ncol = l)
        
        # predefine vec
        vec <- c(1:nrow(h))
        
        for (i in seq(1:nrow(h))){
            vec_temp <- vec[-i]
            numerator_mat[i,  vec_temp] <- p[i, vec_temp]/ (n - m_loc[i] - 
                                                                m_loc[vec_temp] + m_ij[i,  vec_temp])
        }
        
        numerator <- sum(numerator_mat, na.rm = TRUE)
        
        # denominator-------------------------------------------------------------------
        denominator_mat <- matrix(rep(0, l*l), ncol = l)
        
        nullmat <- matrix(rep(1, n*n), ncol=n)
        diag(nullmat) <- 0
        q <- h %*% (nullmat %*% t(h))
        
        for (i in seq(1:nrow(h))){
            vec_temp <- vec[-i]
            denominator_mat[i, vec_temp] <- q[i, vec_temp]/((n - 1) * (n - m_loc[i] - m_loc[vec_temp]) + 
                                            m_loc[i] * m_loc[vec_temp] - m_ij[i, vec_temp])
        }
    
        denominator <- sum(denominator_mat, na.rm = TRUE)
        
        g2_emp <- (numerator / denominator) - 1
        if (verbose == TRUE) {
            if (perm %% 20 == 0) {
                cat("\n", perm, "permutations done")
            } else if (perm == nperm-1) {
                cat("\n", "### permutations finished ###")
            }
            
            if (boot %% 20 == 0) {
                cat("\n", boot, "bootstraps done")
            } else if (boot == nboot-1) {
                cat("\n", "### bootstrapping finished, hell yeah!! ###")
            }
        }
        g2_emp
    }
    
    # g2 point estimate
    g2_emp <- calc_g2(origin)
    
    # permutation of genotypes
    g2_permut <- rep(NA, nperm)
    p_permut <- NA
    
    if (nperm > 0) {
        #setkey(origin, eval(parse(names(origin)[1])))
        perm_genotypes <- function(perm, origin) {
            # origin_perm <- origin[, lapply(.SD, sample)]
            origin_perm <- t(apply(origin, 1, sample)) # to optimize
            g2 <- calc_g2(origin_perm, perm = perm)
            
        }
        if (nperm == 1) nperm <- 2
        g2_permut <- c(g2_emp, sapply(1:(nperm-1), perm_genotypes, origin = origin))
        p_permut <- sum(g2_permut >= g2_emp) / nperm
        perm <- 1
        
    }
    g2_boot <- rep(NA, nboot)
    g2_se <- NA
    CI_boot <- c(NA,NA)
    
    if (nboot > 0) {
        
        boot_genotypes <- function(boot, origin) {
            # bootstrap over individuals in columns
            origin_boot <- origin[, sample(1:ncol(origin), replace = TRUE)]
            g2 <- calc_g2(origin_boot, boot = boot)
        }
        if (nboot == 1) nboot <- 2
        g2_boot <- c(g2_emp, sapply(1:(nboot-1), boot_genotypes, origin = origin))
        g2_se <- stats::sd(g2_boot)
        CI_boot <- stats::quantile(g2_boot, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
    }
    
    res <- list(call=match.call(),
                g2 = g2_emp, p_val = p_permut, g2_permut = g2_permut,
                g2_boot = g2_boot, CI_boot = CI_boot, g2_se = g2_se,
                nobs = nrow(genotypes), nloc = ncol(genotypes))
    
    class(res) <- "inbreed"
    return(res)
    
}
