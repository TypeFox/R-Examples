##' Create genotypes from haplotypes
##'
##' Given a matrix of haplotypes (with one allele per column), this
##'     function creates a matrix of genotypes by combining the
##'     alleles from odd and even columns into genotypes.
##' @title Create genotypes from haplotypes
##' @param haplotypes Matrix of haplotypes (rows: individuals,
##'     columns: alleles)
##' @return Matrix of genotypes (rows: individuals, columns genotypes)
##' @author Karim Oualkacha
##' @keywords internal
Get.G <- function(haplotypes)
{
    m <- dim(haplotypes)[2]/2
    G  <-  NULL
    for(i in 1:m) {
        G <- cbind(G,
                   apply(cbind(haplotypes[, ((2 * i - 1):(2 * i))]),
                         1,
                         ff)
                   )
    }
    return(G)
}


ff <- function(yy) sum(yy[1:2] == 1)
