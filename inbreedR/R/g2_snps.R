#' Estimating g2 from larger datasets, such as SNPs
#'
#' @param genotypes \code{data.frame} with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote), 1 (heterozygote) and NA (missing)
#' @param nperm number or permutations for to estimate a p-value
#' @param nboot number of bootstraps to estimate a confidence interval
#' @param CI confidence interval (default to 0.95)
#' @param parallel Default is FALSE. If TRUE, bootstrapping and permutation tests are parallelized 
#' @param ncores Specify number of cores to use for parallelization. By default,
#'        all available cores are used.
#' @param verbose If FALSE, nothing will be printed to show the status of bootstraps and permutations.
#' 
#'
#' @details Calculates g2 from SNP datasets. Use convert_raw to convert raw genotypes (with 2 columns per locus) into
#'          the required format
#'
#' @return
#' g2_snps returns an object of class "inbreed".
#' The functions `print` and `plot` are used to print a summary and to plot the distribution of bootstrapped g2 values and CI.
#' 
#' An `inbreed` object from \code{g2_snps} is a list containing the following components:
#' \item{call}{function call.}
#' \item{g2}{g2 value}
#' \item{p_val}{p value from permutation test}
#' \item{g2_permut}{g2 values from permuted genotypes}
#' \item{g2_boot}{g2 values from bootstrap samples}
#' \item{CI_boot}{confidence interval from bootstrap distribution}
#' \item{se_boot}{standard error of g2 from bootstraps}
#' \item{nobs}{number of observations}
#' \item{nloc}{number of markers}
#'
#' @references
#' Hoffman, J.I., Simpson, F., David, P., Rijks, J.M., Kuiken, T., Thorne, M.A.S., Lacey, R.C. & Dasmahapatra, K.K. (2014) High-throughput sequencing reveals inbreeding depression in a natural population.
#' Proceedings of the National Academy of Sciences of the United States of America, 111: 3775-3780. Doi: 10.1073/pnas.1318945111
#'
#' @author Martin A. Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Mareike Esser (messer@@techfak.uni-bielefeld.de)
#'
#' @examples
#' # load SNP genotypes in 0 (homozygous), 1 (heterozygous), NA (missing) format.
#' # low number of bootstraps and permutations for computational reasons.
#' data(mouse_snps)
#' (g2_mouse <- g2_snps(mouse_snps, nperm = 10, nboot = 10, CI = 0.95))
#' 
#' # parallelized version for more bootstraps or permutations
#' \dontrun{
#' (g2_mouse <- g2_snps(mouse_snps, nperm = 1000, nboot = 1000, 
#'                      CI = 0.95, parallel = TRUE, ncores = 4))
#' }
#' 
#' @import data.table
#' @export

g2_snps <- function(genotypes, nperm = 0, nboot = 0, CI = 0.95, parallel = FALSE, ncores = NULL, verbose = TRUE) { 
    
    # transpose for congruency with formulae in paper
    origin <- data.table::as.data.table(t(genotypes))
    
    rm(genotypes)
    
    # replace with -1
    for (cols in names(origin)) {
        data.table::set(origin, i=which((origin[[cols]] != 1L & origin[[cols]] != 0L) | 
                is.na(origin[[cols]] != 0L)), j=cols, value= -1L)
    }
    
    n <- ncol(origin) # number of individuals
    l <- nrow(origin) # number of loci
    
    
    calc_g2 <- function(origin, perm = 1, boot = 1) {
        
        # H matrix with 0 for -1
        h <- data.table::copy(origin)
        
        for (cols in seq_along(names(h))) {
            data.table::set(h, i=which((h[[cols]] != 1L)  & (h[[cols]] != 0L)), 
                j=cols, value= 0L)
        }
        
        # precalculating sums
        rowsum_h <- rowSums(h, na.rm = TRUE)
        colsum_h <- colSums(h, na.rm = TRUE)
        sum_rowsum_squared <- sum(rowsum_h^2)
        sum_colsum_squared <- sum(colsum_h^2)
        h_sum <- sum(colsum_h, na.rm = TRUE)
        
        rm(h)
        
        # define matrix with 1 for missing and 0 for all others
        m <- data.table::copy(origin)
        
        for (cols in seq_along(names(m))) {
            data.table::set(m, i=which((m[[cols]] == 1L)), j=cols, value= 0L)
            data.table::set(m, i=which((m[[cols]] == -1L)), j=cols, value= 1L)
        }
        
        # vector with rowsums for missing data matrix
        m_loc_temp <- rowSums(m, na.rm = TRUE)
        m_loc <- m_loc_temp / n
        
        
        # numerator
        numer <- (n-1) * (sum_colsum_squared - h_sum) /
            (h_sum^2 - sum_rowsum_squared - sum_colsum_squared + h_sum)
        
        # overwrite m
        for (cols in seq_along(names(m))) {
            data.table::set(m,
                j= cols, 
                value= (rowsum_h * m[[cols]]) / (1 - m_loc))
        }
        
        # delete infinity
        for (j in seq_along(names(m))) {
            data.table::set(m, 
                i = which(is.infinite(m[[j]])), 
                j = j,
                value = NA)
        }
        
        M_ind <- colSums(m, na.rm = TRUE)^2 - colSums(m^2, na.rm = TRUE)
        
        # delete missmat
        rm(m)
        
        M_hat <- (1/(n)) * sum(M_ind, na.rm = TRUE)
        X_temp <- rowsum_h * m_loc / (1-m_loc)
        
        a_hat_temp <- (sum(X_temp^2, na.rm = TRUE) - (sum(X_temp, na.rm = TRUE))^2)
        a_hat <- (M_hat + a_hat_temp) / (h_sum^2 - sum_rowsum_squared)
        g2_emp <- numer / (1 + a_hat) - 1
        
        if (verbose == TRUE){
            if (perm %% 5 == 0) {
                cat("\n", perm, "permutations done")
            } else if (perm == nperm-1) {
                cat("\n", "### permutations finished ###")
            }
            
            if (boot %% 5 == 0) {
                cat("\n", boot, "bootstraps done")
            } else if (boot == nboot-1) {
                cat("\n", "### bootstrapping finished ###")
            }
        }
        g2_emp
    }
    # g2 point estimate
    g2_emp <- calc_g2(origin)
    
    # permutation of genotypes
    g2_permut <- rep(NA, nperm)
    p_permut <- NA
    
    # permutation function
    perm_genotypes <- function(perm, origin) {
        # columnwise permutation
        origin_perm <- origin[, lapply(.SD, sample)]
        # origin_perm <- lapply(origin, sample)
        g2 <- calc_g2(origin_perm, perm = perm)
    }
    
    # 1 gets substracted due to the adding of the empirical value
    if (nperm == 1) nperm <- 2 
    
    if (nperm > 0 & parallel == FALSE) {
        
        g2_permut <- sapply(1:(nperm-1), perm_genotypes, origin)
        p_permut <- sum(c(g2_emp, g2_permut) >= g2_emp) / nperm
        perm <- 1
        
    }
    
    if (nperm > 0 & parallel == TRUE) {
        if (is.null(ncores)) {
            ncores <- parallel::detectCores()
            warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
        }
        
        cl <- parallel::makeCluster(ncores)
        g2_permut <- parallel::parSapply(cl, 1:(nperm-1), perm_genotypes, origin)
        parallel::stopCluster(cl)
        p_permut <- sum(c(g2_emp, g2_permut) >= g2_emp) / nperm
        perm <- 1
    }
    
    # bootstrap
    g2_boot <- rep(NA, nboot)
    g2_se <- NA
    CI_boot <- c(NA,NA)
    
    
    # bootstap function
    boot_genotypes <- function(boot, origin) {
        origin_boot <- origin[, sample(1:ncol(origin), replace = TRUE), with = FALSE]
        g2 <- calc_g2(origin_boot, boot = boot)
    }
    
    if (nboot == 1) nboot <- 2
    
    if (nboot > 0 & parallel == FALSE) {
        g2_boot <- c(g2_emp, sapply(1:(nboot-1), boot_genotypes, origin))
        g2_se <- stats::sd(g2_boot)
        CI_boot <- stats::quantile(g2_boot, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
    }
    
    if (nboot > 0 & parallel == TRUE) {
        if (is.null(ncores)) {
            ncores <- parallel::detectCores()
            warning("No core number specified: detectCores() is used to detect the number of \n cores on the local machine")
        }
        
        # start cluster
        cl <- parallel::makeCluster(ncores)
        g2_boot <- c(g2_emp, parallel::parSapply(cl, 1:(nboot-1), boot_genotypes, origin))
        parallel::stopCluster(cl)
        g2_se <- stats::sd(g2_boot)
        CI_boot <- stats::quantile(g2_boot, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
        
        #             R.boot <- unname(parallel::parApply(cl, Ysim, 2, R.pe, groups = groups))
        #             parallel::stopCluster(cl)
        
    }
    
    res <- list(call=match.call(),
        g2 = g2_emp, p_val = p_permut, g2_permut = g2_permut,
        g2_boot = g2_boot, CI_boot = CI_boot, g2_se = g2_se,
        nobs = ncol(origin), nloc = nrow(origin))
    class(res) <- "inbreed"
    return(res)
}