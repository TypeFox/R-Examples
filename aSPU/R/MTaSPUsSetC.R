#' gene-Multitrait Sum of Powered Score (MTSPUsSet) tests and adaptive MTSPUsSet (MTaSPUsSet) test for multi trait - SNP set association with GWAS summary statistics. (C coded version)
#'
#' It gives p-values of the MTSPUsSet tests and MTaSPUsSet test with GWAS summary statistics.
#'
#' @param Zs Z-scores for each SNPs. It could be P-values if the Ps option is TRUE. 
#'
#' @param corSNP Correlation matirx of the SNPs to be tested; estimated from a
#' reference panel (based on the same set of the reference alleles as
#' used in calculating Z-scores).
#'
#' @param corPhe Correlation matirx of phenotypes to be tested; Estimated from Z-scores.
#' 
#' @param pow SNP specific power(gamma values) used in MTSPUsSet test.
#'
#' @param pow2 GENE specific power(gamma values) used in MTSPUsSet test.
#'
#' @param n.perm number of permutations or bootstraps.
#'
#' @param Ps TRUE if input is p-value, FALSE if input is Z-scores. The default is FALSE.
#'
#' @return A vector object, MTSPUsSet test P values and MTaSPUsSet P value.
#'
#' @author Il-Youp Kwak and Wei Pan
#'
#' @references
#'
#' Il-Youp Kwak, Wei Pan (2016)
#'Gene- and pathway-based association tests for multiple
#'      traits with GWAS summary statistics
#'
#' @examples
#'
#' data(SAMD11)
#' attach(SAMD11)
#' ## example analysis using aSPUM test.
#' (outFZC <- MTaSPUsSetC(ZsF, corSNP=corSNPF, corPhe = corPheF,
#'       pow=c(1,2,4,8),  pow2 = c(1,2,4,8), n.perm=10, Ps=FALSE))
#'
#' 
#' @seealso \code{\link{MTaSPUsSet}} 

MTaSPUsSetC <- function(Zs, corSNP, corPhe, pow=c(1,2,4,8),
                     pow2 = c(1,2,4,8),
                     n.perm=5000,
                     Ps = FALSE) {
    if(!is.matrix(Zs))
        Zs <- t(as.matrix(Zs))
    nsnp <- dim(Zs)[1]
    nphe <- dim(Zs)[2]

    V <- corSNP
    U <- corPhe

    e.U <- eigen(U)
    e.U$values[ e.U$values < 0 ] = 0
    A <- e.U$vectors %*% diag(sqrt(e.U$values)) %*% t(e.U$vectors)
    e.V <- eigen(V)
    e.V$values[ e.V$values < 0 ] = 0
    B <- e.V$vectors %*% diag(sqrt(e.V$values)) %*% t(e.V$vectors)

    Zs <- as.matrix(Zs)
    if(Ps == TRUE)
        Zs <- qnorm(1 - Zs/2)

    Ts <- rep(0, length(pow)*length(pow2))

    for(p1 in 1:length(pow)) {
        for(p2 in 1:length(pow2)) {

            if(pow[p1] < Inf) {
                Ts[p2 + (p1-1) * length(pow2)] = sum(apply(Zs, 2, function(x) sum(x^pow[p1]) )^pow2[p2])
            } else {
                Ts[p2 + (p1-1) * length(pow2)] = sum(apply(Zs, 2, max )^pow2[p2])
            }
        }
    }

    ## Permutations:

    pow[pow==Inf] = 0 # pass 0 as infitiy
    if(Ps == TRUE)
        Pval = 1 else Pval = 2

    T0sC <- calcT0simM2(as.matrix(A), as.matrix(B), as.matrix(pow), as.matrix(pow2), n.perm, Pval,t(as.matrix(Ts)) )

                                        #    minP0s <- T0sC$mPs/n.perm
                                        #    minP =  sum( min(T0sC$Pvs) > minP0s )/n.perm

    pvs=c(T0sC$Pvs, T0sC$minp/n.perm)
    ##    pvs=c(T0sC$Pvs, minp)

    nmvec <- NULL;
    for(ii in pow2) {
        for(jj in pow) {
            nmvec <- c(nmvec, paste("MTSPUs",jj,",",ii,sep=""))
        }
    }

    nmvec <- c(nmvec, "MTaSPUs")
    names(pvs) <- nmvec

    pvs
}
