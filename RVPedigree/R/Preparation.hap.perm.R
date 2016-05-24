##' Prepare for the haplotype permutation
##'
##' This functions prepares for running the \code{\link{perm.hap}}
##' function. The output of this function can be passed to that one.
##' @param pedigree data frame with pedigree information in linkage
##' format (four columns: family ID, individual ID, father ID, mother
##' ID). For example the first for columns of a Plink \code{.ped} file
##' or a Plink \code{.bim} file.
##' @param generation.id a vector of length(sample size) which
##'     indicates if the subject is founder (\code{generation.id=0}),
##'     a child from first generation (\code{generation.id=1}), a
##'     child from second generation (\code{generation.id=2}), etc.
##'     This vector can be calculated by the
##'     \code{kinship2::kindepth()} function.
##' @param indices.haplotypes1
##' @param indices.haplotypes2
##' @return a list of vectors with indices of various parts of the
##' pedigree.
##' @seealso \code{\link{read.pedigree}}
##' @author Karim Oualkacha
##' @author M'Hamed Lajmi Lakhal-Chaieb
##' @keywords internal
Preparation.hap.perm <- function(pedigree,
                                 generation.id,
                                 indices.haplotypes1,
                                 indices.haplotypes2){
    indices   <- 1:nrow(pedigree)
    fam.id    <- pedigree[, 1]
    ind.id    <- pedigree[, 2]
    father.id <- pedigree[, 3]
    mother.id <- pedigree[, 4]

    is.founder       <- (father.id==0)
    indices.founders <- indices[is.founder]
    n.founders       <- sum(is.founder)
    n.hap.founders   <- 2 * n.founders

    Result1     <- NULL
    Result      <- NULL
    indices.gen <- NULL

    nbre.generations <- length(unique(generation.id)) - 1
    n.gen            <- rep(-999, nbre.generations)

    for (i in 1:nbre.generations)
    {
        # A list with TRUE for all people in generation i (and not founder)
        is.gen           <- ((generation.id == i) & (father.id != 0))
        # A vector with the indices of all people in generation i
        indices.gen[[i]] <- indices[is.gen]
        # The number of people in generation i
        n.gen[i]         <- sum(is.gen)
        # The indices of the parents of the people in generation i
        indices.father   <- sapply(indices.gen[[i]],
                                   findParentIndices,
                                   indices=indices,
                                   fam.id=fam.id,
                                   ind.id=ind.id,
                                   parent.id=father.id)
        indices.mother   <- sapply(indices.gen[[i]],
                                   findParentIndices,
                                   indices=indices,
                                   fam.id=fam.id,
                                   ind.id=ind.id,
                                   parent.id=mother.id)
        # An array with the indices of the parents as columns
        indices.fam      <- cbind(indices.father, indices.mother)
        Result1[[i]]     <- indices.fam  # indices of families: indices.fam1 = sss$indices.fam1; indices.fam2=sss$indices.fam2;
    }

    Result$n.founders          <- n.founders
    Result$indices.founders    <- indices.founders
    Result$indices.fam         <- Result1
    Result$n.gen               <- n.gen  # n.gen1=sss$n.gen1;n.gen2=sss$n.gen2;
    Result$indices.gen         <- indices.gen
    Result$n.hap.founders      <- 2 * n.founders
    Result$indices.haplotypes1 <- indices.haplotypes1
    Result$indices.haplotypes2 <- indices.haplotypes2


    return(Result)
}
