##' Compute the p-value for the VC-C2 method
##'
##' @inheritParams RVPedigree
##' @param P a matrix which can be obtained as output of the function
##' \code{\link{Preparation.VCC}}.
##' @param G the genotype matrix
##' @param n Integer, the number of samples
##' @param pedigree a pedigree as output by
##' \code{\link{read.pedigree}}. This is a data frame consisting of
##' four columns (family ID, individual ID, father ID and mother ID)
##' as use in the traditional linkage format and e.g. Plink files.
##' @param haplotypes a matrix of the haplotypes of the individuals
##' @param generation.id a vector of length(sample size) which
##'     indicates if the subject is founder (\code{generation.id=0}),
##'     a child from first generation (\code{generation.id=1}), a
##'     child from second generation (\code{generation.id=2}), etc.
##'     This vector can be calculated by the
##'     \code{kinship2::kindepth()} function.
##' @param Ncores (integer) Number of processor (CPU) cores to be used
##'     in parallel when doing the permutations to determine the
##'     p-value (default: 1).
##' @return a list with the score and p-value for the VC-C2
##'     association test on the given region. The list contains the
##'     elements \code{score} and \code{p.value}. This is a
##'     permutation-based p-value.
##' @author Karim Oualkacha
##' @author M'Hamed Lajmi Lakhal-Chaieb
##' @seealso \code{\link{pvalue.VCC1}}, \code{\link{pvalue.VCC3}}
##' @keywords internal
##' @import foreach doParallel
pvalue.VCC2 <- function(P, G, W, Nperm, n, pedigree, haplotypes,
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

    res.perm <- foreach(b=1:Nperm, .combine='c', .verbose=FALSE) %dopar% {
        haplotypes.perm <- perm.hap(haplotypes,
                                    preplist=prep.hap,
                                    generation.id=generation.id)
        G.perm          <- Get.G(haplotypes.perm)
        sum(diag(t(G.perm) %*% P %*% G.perm %*% W))
    }

    p.value <- mean(res.perm>=score)
    return(list(score=score,
                p.value=p.value))
}
