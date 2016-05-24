#' Expected r2 between standardized multilocus heterozygosity (h) and inbreeding level (f)
#' 
#' 
#'
#' @param genotypes \code{data.frame} with individuals in rows and loci in columns,
#'        containing genotypes coded as 0 (homozygote), 1 (heterozygote) and \code{NA} (missing)
#' @param type specifies g2 formula to take. Type "snps" for large datasets and "msats" for smaller datasets.
#' @param nboot number of bootstraps over individuals to estimate a confidence interval
#'        around r2(h, f)
#' @param parallel Default is FALSE. If TRUE, bootstrapping and permutation tests are parallelized 
#' @param ncores Specify number of cores to use for parallelization. By default,
#'        all available cores but one are used.
#' @param CI confidence interval (default to 0.95)
#' 
#' @return 
#' \item{call}{function call.}
#' \item{r2_hf_full}{expected r2 between inbreeding and sMLH for the full dataset}
#' \item{r2_hf_boot}{expected r2 values from bootstrapping over individuals}
#' \item{CI_boot}{confidence interval around the expected r2}
#' \item{nobs}{number of observations}
#' \item{nloc}{number of markers}
#' 
#' 
#' @references
#' Slate, J., David, P., Dodds, K. G., Veenvliet, B. A., Glass, B. C., Broad, T. E., & McEwan, J. C. (2004). 
#' Understanding the relationship between the inbreeding coefficient 
#' and multilocus heterozygosity: theoretical expectations and empirical data. Heredity, 93(3), 255-265.
#' 
#' Szulkin, M., Bierne, N., & David, P. (2010). HETEROZYGOSITY-FITNESS CORRELATIONS: A TIME FOR REAPPRAISAL. 
#' Evolution, 64(5), 1202-1217.
#' 
#' @author Martin A. Stoffel (martin.adam.stoffel@@gmail.com) 
#'        
#' @examples
#' data(mouse_msats)
#' genotypes <- convert_raw(mouse_msats)
#' (out <- r2_hf(genotypes, nboot = 100, type = "msats", parallel = FALSE))
#' plot(out)
#' @export
#'
#'

r2_hf <- function(genotypes, type = c("msats", "snps"), nboot = NULL, 
                  parallel = FALSE, ncores = NULL, CI = 0.95) {
    
        
#     if (!(steps > 1) | (steps > ncol(genotypes))) {
#         stop("steps have to be at least two and smaller or equal than the number of markers used")
#     }
    gtypes <- as.matrix(genotypes)
    
    # check g2 function argument
    if (length(type) == 2){
        type <- "msats"
    } else if (!((type == "msats")|(type == "snps"))){
        stop("type argument needs to be msats or snps")
    } 
    
    # define g2 function
    if (type == "msats") {
        g2_fun <- g2_microsats
    } else if (type == "snps") {
        g2_fun <- g2_snps
    }
    
    # define calculation of expected r2
    calc_r2 <- function(gtypes) {
            g2 <- g2_fun(gtypes)[["g2"]]
            # according to the miller paper, negative g2s are set to r2 = 0.
            if (g2 < 0) return(r2 <- 0)
            var_sh <- stats::var(sMLH(gtypes))
            r2 <- g2/var_sh
            if (r2 > 1) r2 <- 1
            r2
    }
    
    # calculate r2 for full data
    r2_hf_full <- calc_r2(gtypes)
    
    # initialise r2_hf_boot
    r2_hf_boot <- NA
    CI_boot <- NA
    # bootstrapping over r2?
    if (!is.null(nboot)) {
        if (nboot < 2) stop("specify at least 2 bootstraps with nboot to 
                            estimate a confidence interval")
        # initialise
        r2_hf_boot <- matrix(nrow = nboot)
        
        # bootstrap function
        calc_r2_hf_boot <- function(boot, gtypes) {
            inds <- sample(1:nrow(gtypes), replace = TRUE)
            out <- calc_r2(gtypes[inds, ])
            # notifications
            if (boot %% 20 == 0) cat("\n", boot, "bootstraps over individuals done")
            if (boot == nboot) cat("\n", "### bootstrapping over individuals finished! ###")
            # result
            return(out)
        }
        
        # parallel bootstrapping ?
        if (parallel == FALSE) r2_hf_boot <- sapply(1:nboot, calc_r2_hf_boot, gtypes)
        
        if (parallel == TRUE) {
            if (is.null(ncores)) {
                ncores <- parallel::detectCores()-1
                warning("No core number specified: detectCores() is used to detect the number 
                        of \n cores on the local machine")
            }
            # run cluster
            cl <- parallel::makeCluster(ncores)
            parallel::clusterExport(cl, c("g2_microsats", "g2_snps", "sMLH"), envir = .GlobalEnv)
            r2_hf_boot <- parallel::parSapply(cl, 1:nboot, calc_r2_hf_boot, gtypes)
            parallel::stopCluster(cl)
        }
        r2_hf_boot <- c(r2_hf_full, r2_hf_boot)
        CI_boot <- stats::quantile(r2_hf_boot, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
    }
    
 
    
    res <- list(call = match.call(),
                r2_hf_full = r2_hf_full,
                r2_hf_boot = r2_hf_boot,
                CI_boot = CI_boot,
                nobs = nrow(genotypes), 
                nloc = ncol(genotypes))
    
    class(res) <- "inbreed"
    return(res)
    
}