#' GRAMMAR test for association in samples with genetic structure
#'
#' Fast approximate test for association between a trait and genetic polymorphisms,
#' in samples with genetic sub-structure (e.g. relatives). The function implements
#' several varieties of GRAMMAR ('gamma','gc', and 'raw').
#'
#' With 'raw' argument,
#' the original GRAMMAR (Aulchenko et al., 2007) is implemented. This method
#' is conservative and generates biased estimates of regression coefficients.
#'
#' With 'gc' argument, the GRAMMAR-GC (Amin et al., 2007) is implemented.
#' This method solves the conservativity of the test, but the Genomic Control (GC)
#' lambda is by definition "1" and can not serve as an indicator of goodness of
#' the model; also, the estimates of regression coefficients are biased (the
#' same as in 'raw' GRAMMAR).
#'
#' GRAMMAR-Gamma (default 'gamma' argument) solves these problems, producing
#' a correct distribution of the test statistic, interpretable value of GC Lambda,
#' and unbiased estimates of the regression coefficients. All together, the
#' default 'gamma' method is recommended for use.
#'
#' @param polyObject object returned by \code{\link{polygenic}} function
#' @param data object of \code{\link{gwaa.data-class}}
#' @param method to be used, one of 'gamma','gc', or 'raw'
#' @param propPs proportion of non-corrected P-values used to estimate
#' the inflation factor Lambda, passed directly to the \code{\link{estlambda}}
#' @param ... arguments passed to the function used for computations,
#' (\code{\link{qtscore}})
#'
#' @return Object of scan.gwaa-class
#'
#' @seealso \code{\link{polygenic}}, \code{\link{mmscore}}, \code{\link{qtscore}}
#'
#' @references
#'
#' GRAMMAR-Raw:
#' Aulchenko YS, de Koning DJ, Haley C.
#' Genomewide rapid association using mixed model and regression: a fast and
#' simple method for genomewide pedigree-based quantitative trait loci
#' association analysis. Genetics. 2007 Sep;177(1):577-85.
#'
#' GRAMMAR-GC:
#' Amin N, van Duijn CM, Aulchenko YS.
#' A genomic background based method for association analysis in related individuals.
#' PLoS One. 2007 Dec 5;2(12):e1274.
#'
#' GRAMMAR-Gamma:
#' Svischeva G, Axenovich TI, Belonogova NM, van Duijn CM, Aulchenko YS.
#' Rapid variance components-based method for whole-genome association analysis.
#' Nature Genetics. 2012 44:1166-1170. doi:10.1038/ng.2410
#'
#' @examples
#' # Using clean ge03d2 data
#' require(GenABEL.data)
#' data(ge03d2.clean)
#' # take only a small piece for speed
#' ge03d2.clean <- ge03d2.clean[1:200,]
#' # estimate genomic kinship
#' gkin <- ibs(ge03d2.clean[,sample(autosomal(ge03d2.clean),1000)], w="freq")
#' # perform polygenic analysis
#' h2ht <- polygenic(height ~ sex + age, kin=gkin, ge03d2.clean)
#' h2ht$est
#' # compute mmscore stats
#' mm <- mmscore(h2ht, data=ge03d2.clean)
#' # compute grammar-gc
#' grGc <- grammar(h2ht, data=ge03d2.clean, method="gc")
#' # compute grammar-gamma
#' grGamma <- grammar(h2ht, data=ge03d2.clean, method="gamma")
#' # compare lambdas
#' lambda(mm)
#' estlambda(mm[,"chi2.1df"])
#' lambda(grGamma)
#' estlambda(grGamma[,"chi2.1df"])
#' lambda(grGc)
#' estlambda(grGc[,"chi2.1df"])
#' # compare top results
#' summary(mm)
#' summary(grGamma)
#' summary(grGc)
#'
#' @author Gulnara Svischeva, Yurii Aulchenko
#'
#' @keywords htest
#'
"grammar" <-
    function(polyObject, data, method=c("gamma","gc","raw"), propPs=1.0, ... )
{
    method <- match.arg(method)
    if (method == "gamma") {
        out <- qtscore(polyObject$pgresidualY, data=data, clambda=TRUE, ... )

        ## correct test and beta values
        out@results[,"chi2.1df"] <-
            out@results[, "chi2.1df"]/polyObject$grammarGamma$Test
        out@results[,"effB"] <- out@results[, "effB"]/polyObject$grammarGamma$Beta

        ## re-compute standard errors (c2 = (b/se)^2; se = b/sqrt(c2))
        out@results[,"se_effB"] <-
            abs(out@results[, "effB"])/sqrt(out@results[, "chi2.1df"])
        ## recompute p-values
        out@results[,"P1df"] <-
            pchisq(out@results[,"chi2.1df"], df=1, lower.tail=FALSE)
        ## recompute Lambda
        out@lambda <- estlambda(out[, "chi2.1df"],
                                plot=FALSE,
                                proportion=propPs)

        if (out@lambda$estimate <= 1) {
            warning("estimate of Lambda < 1, constraining to 1")
            out@lambda$estimate <- 1.0
            out@lambda$se <- NA
        }
    } else if (method == "gc") {
        out <- qtscore(polyObject$pgres, data=data,
                       clambda=FALSE, propPs=propPs, ... )
    } else if (method == "raw") {
        out <- qtscore(polyObject$pgres, data=data,
                       clambda=TRUE, propPs=propPs, ... )
    } else {
        stop("method should be one of 'gamma','gc', or 'raw'")
    }

    ## set uncorrected stats to NA to avoid confusion
    out@results[,"effAB"] <- out@results[,"effBB"] <- out@results[,"chi2.2df"] <-
        out@results[,"P2df"] <- NA
    out@call <- match.call()
    ## return results
    return(out);
}
