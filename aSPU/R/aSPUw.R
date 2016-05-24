#' Inverse variance weighted Sum of Powered Score tests (SPUw) and adaptive SPUw (aSPUw) test for single trait - SNP set association.
#'
#' It gives the p-values of the SPUw tests and aSPUw test based
#' on the permutations of the residuals or simulations from the null distripution.
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
#' @export
#' @return A list object, Ts : Test Statistics for the SPUw and aSPUw test.
#'         pvs : p-values for the SPUw and aSPUw test.
#'
#' @author Il-Youp Kwak, Junghi Kim and Wei Pan
#'
#' @references
#' Junghi Kim, Jeffrey R Wozniak, Bryon A Mueller, Xiaotong Shen and Wei Pan (2014)
#' Comparison of statistical tests for group differences in brain functional networks,
#' Neuroimage, 1;101:681-694
#'
#' @examples
#'
#' data(exdat)
#' out <- aSPUw(exdat$Y, exdat$X, pow = c(1:8, Inf), n.perm = 1000)
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
#' @seealso \code{\link{aSPU}}


aSPUw <- function(Y, X, cov=NULL, resample = c("perm", "sim", "boot"), model=c("gaussian", "binomial"), pow = c(1:8, Inf), n.perm = 1000) {

    model <- match.arg(model)
    resample <- match.arg(resample)
    userank = T

    if(resample == "sim") {
        aSPUwsim(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)
    } else {
        if (resample == "boot") {
            aSPUwboot(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model)
	} else {
        aSPUwpermC(Y = Y, X = X, cov = cov, pow = pow, n.perm = n.perm, model = model, userank = userank)
    	}
    }
}
