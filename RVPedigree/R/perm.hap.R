##' Permute haplotypes
##'
##' @param haplotypes a matrix of haplotypes to be permuted
##' @param preplist a list as output by \code{\link{Preparation.hap.perm}}
##' @param generation.id a vector of length(sample size) which
##'     indicates if the subject is founder (\code{generation.id=0}),
##'     a child from first generation (\code{generation.id=1}), a
##'     child from second generation (\code{generation.id=2}), etc.
##'     This vector can be calculated by the
##'     \code{kinship2::kindepth()} function.
##' @return a matrix of permuted haplotypes
##' @author Karim Oualkacha
##' @author M'Hamed Lajmi Lakhal-Chaieb
##' @keywords internal
perm.hap <- function(haplotypes, preplist, generation.id)
{
    ## haplotypes is a matrix of subjects haplotypes

    ## NB: we need only founders haplotypes. In the running example,
    ## since it is a simualted data, we know haplotypes for all
    ## individuals, thus, we used them all.

    haplotypes.perm <- matrix(-999,
                              ncol=ncol(haplotypes),
                              nrow=nrow(haplotypes))

    haplotypes.founders <- rbind(haplotypes[preplist$indices.founders,
                                            preplist$indices.haplotypes1],
                                 haplotypes[preplist$indices.founders,
                                            preplist$indices.haplotypes2])

    haplotypes.founders.perm <-
        haplotypes.founders[sample(1:preplist$n.hap.founders), ]

    haplotypes.perm[preplist$indices.founders,
                    preplist$indices.haplotypes1] <-
        haplotypes.founders.perm[(1:preplist$n.founders), ]

    haplotypes.perm[preplist$indices.founders,
                    preplist$indices.haplotypes2] <-
        haplotypes.founders.perm[((preplist$n.founders + 1):preplist$n.hap.founders), ]


    nbre.generations <- length(unique(generation.id)) - 1

    for (i in 1:nbre.generations){
        haplotypes.perm[preplist$indices.gen[[i]],
                        preplist$indices.haplotypes1] <-
            t(sapply((1:preplist$n.gen[i]),
                     meiose,
                     indices.fam=preplist$indices.fam[[i]],
                     indices.haplotypes=preplist$indices.haplotypes1,
                     haplotypes.perm=haplotypes.perm))

        haplotypes.perm[preplist$indices.gen[[i]],
                        preplist$indices.haplotypes2] <-
            t(sapply((1:preplist$n.gen[i]),
                     meiose,
                     indices.fam=preplist$indices.fam[[i]],
                     indices.haplotypes=preplist$indices.haplotypes2,
                     haplotypes.perm=haplotypes.perm))
    }

    return(haplotypes.perm)
}
