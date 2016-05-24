##' Generate permutations of a phenotype vector
##' 
##' Given a vector, generate n.perm samples and return a matrix with each
##' permutation in each column.
##' 
##' 
##' @param pheno a vector to be permuted
##' @param n.perm the number of times to permute
##' @return a matrix with dimensions length(pheno) x n.perm.
##' @author Chris Wallace <chris.wallace at cimr.cam.ac.uk>
##' @export
##' @keywords manip
##' @examples
##' 
##' y <- rbinom(50,2,0.3)
##' genperms(y,4)
##' 
genperms <-
function(pheno,n.perm=0) {
  pheno.perm <- matrix(as.integer(0),length(pheno),n.perm)
  for(j in 1:n.perm) {
    pheno.perm[,j] <- sample(pheno)
  }
  return(pheno.perm)
}
