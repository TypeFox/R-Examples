##' Compute p-values for the VC-C3 method
##'
##' @inheritParams pvalue.VCC2
##' @return the outcomes of the VC-C3 association test on the given
##'     region, which is a list with the following elements:
##'     \itemize{
##'     \item \code{score}: The score test
##'     \item \code{p.value.VCC2}: The p-value as output by the VC-C2
##'     method (\code{\link{pvalue.VCC2}}).
##'     \item \code{p.value.VCC3}: This is the exact p-vlaue for small
##'     level size. The approach to obtain this second p-value
##'     approximates the distribution of the test statistic under the
##'     null model by a non-central chi-square using matching moments
##'     of the first, the second and the fourth moments. The
##'     estimation of the three moments are approximated empirically
##'     using the permutations that were used to calculate the VC-C2
##'     p-value.
##'     }
##' @author Karim Oualkacha
##' @author M'Hamed Lajmi Lakhal-Chaieb
##' @seealso \code{\link{pvalue.VCC1}}, \code{\link{pvalue.VCC2}}
##' @keywords internal
pvalue.VCC3 <- function(P, G, W, Nperm, n, pedigree, haplotypes,
                        generation.id, Ncores=1)
{
    ## Set up parallel environment for computing the permutations. By
    ## default Ncores=1, so no problems are expected.
    check_Ncores(Ncores)
    registerDoParallel(Ncores)

    nrHaplotypes <- ncol(haplotypes)
    indices.haplotypes1 <- seq(1, nrHaplotypes, 2)
    indices.haplotypes2 <- seq(2, nrHaplotypes, 2)

    prep.hap <- Preparation.hap.perm(pedigree=pedigree,
                                     generation.id=generation.id,
                                     indices.haplotypes1=indices.haplotypes1,
                                     indices.haplotypes2=indices.haplotypes2
                                     )

    score <- sum(diag(t(G) %*% P %*% G %*% W))

    res.perm <- foreach(b=1:Nperm, .combine='c') %dopar% {
        haplotypes.perm <- perm.hap(haplotypes,
                                    preplist=prep.hap,
                                    generation.id=generation.id)
        G.perm <- Get.G(haplotypes.perm)
        sum(diag(t(G.perm) %*% P %*% G.perm %*% W))
    }

    # This is the same p-value as output by pvalue.VCC2()
    p.value.VCC2 <- mean(res.perm >= score)

    m1 <- mean(res.perm)
    m2 <- mean((res.perm - m1)^2)
    m4 <- mean((res.perm - m1)^4)
    g <- (m4 / (m2^2)) - 3
    df <- 12 / g
    stat <- df + (score - m1) * sqrt(2 * df / m2)
    p.value.VCC3 <- 1 - pchisq(stat, df=df)

    return(list(score=score,
                p.value.VCC2=p.value.VCC2,
                p.value.VCC3=p.value.VCC3))
}
