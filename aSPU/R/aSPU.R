#' Sum of Powered Score (SPU) tests and adaptive SPU (aSPU) test for single trait - SNP set association.
#'
#' It gives p-values of the SPU tests and aSPU test.
#'
#' @param Y Response or phenotype data. It can be a disease indicator; =0 for controls, =1 for cases.
#' Or it can be a quantitative trait. A vector with length n (number of observations).
#'
#' @param X Genotype or other data; each row for a subject, and each column
#'     for an SNP (or a predictor). The value of each SNP is the # of the copies
#'     for an allele. A matrix with dimension n by k (n : number of observation, k : number of SNPs (or predictors) ).
#'
#' @param cov Covariates. A matrix with dimension n by p (n :number of observation, p : number of covariates).
#'
#' @param resample Use "perm" for residual permutations, "sim" for simulations from the null distribution, and "boot" for parametric bootstrap.
#'
#' @param model Use "gaussian" for a quantitative trait, and use "binomial" for a binary trait.
#'
#' @param pow power used in SPU test. A vector of the powers.
#'
#' @param n.perm number of permutations or bootstraps.
#'
#' @return A list object, Ts : test statistics for the SPU tests (in the order of the specified pow) and finally for the aSPU test.
#'         pvs : p-values for the SPU and aSPU tests.
#'
#' @author Il-Youp Kwak, Junghi Kim, Yiwei Zhang and Wei Pan
#'
#' @references
#' Wei Pan, Junghi Kim, Yiwei Zhang, Xiaotong Shen and Peng Wei (2014)
#' A powerful and adaptive association test for rare variants,
#' Genetics, 197(4), 1081-95
#'
#' Junghi Kim, Jeffrey R Wozniak, Bryon A Mueller, Xiaotong Shen and Wei Pan (2014) Comparison of statistical tests for group differences in brain functional networks, NeuroImage, 1;101:681-694
#'
#' @examples
#'
#' data(exdat)
#'
#' ## example analysis using aSPU test on exdat data.
#' out <- aSPU(exdat$Y, exdat$X, cov = NULL, resample = "boot",
#'            model = "binomial", pow = c(1:8, Inf), n.perm = 1000)
#'
#' out$Ts
#' # This is a vector of Test Statistics for SPU and aSPU tests.
#' # SPU1 to SPUInf corresponds with the option pow=c(1:8, Inf)
#' # They are SPU test statistics.
#' # The last element aSPU is minimum of them, aSPU statistic.
#'
#' out$pvs
#' # This is a vector of p-values for SPU and aSPU tests.
#' # SPU1 to SPUInf corresponds with the option pow=c(1:8, Inf)
#' # They are p-values for corresponding SPU tests.
#' # The last element is p-value of aSPU test.
#'
#' @seealso \code{\link{aSPUw}}

aSPU <- function(Y, X, cov=NULL, resample = c("perm", "sim", "boot"), model=c("gaussian", "binomial"), pow = c(1:8, Inf), n.perm = 1000) {
    X <- as.matrix(X)
    model <- match.arg(model)
    resample <- match.arg(resample)

    if(resample == "boot") {
        if(n.perm > 10^5) {
            aSPUboot(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)
        } else {
            aSPUboot2(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)
        }
    } else {
        if(resample == "sim") {
            aSPUsimC2(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)            
        } else {
            aSPUpermC2(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)
        }
    }
}
