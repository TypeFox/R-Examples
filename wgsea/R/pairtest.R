##' Generate p values for each SNP for case-control comparisons.
##' 
##' A wrapper for the \code{snpStats} function \code{single.snp.tests}.
##' Generates p values for the association of each SNP with case or control
##' status.
##' 
##' 
##' @param case \code{SnpMatrix} object holding genotypes of case subjects
##' @param control \code{SnpMatrix} object holding genotypes of control
##' subjects
##' @param n.perm number of permutations of case control status required to
##' generate permuted p value vectors.  The default, given by \code{n.perm=0},
##' is not to permute.
##' @param pheno.perm An alternative to specifying \code{n.perm} is to supply a
##' matrix of alternative phenotypes, with each column relating to a different
##' permutation.
##' @param quiet set TRUE to suppress the printing of progress dots
##' @return If \code{n.perm=0}, a vector of p values, one for each SNP (each
##' column in the \code{case} and \code{control} objects.  If \code{n.perm>0},
##' a matrix of p values, each column representing the results of a different
##' permutation.  %% ~Describe the value returned %% If it is a LIST, use %%
##' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
##' 'comp2'} %% ...
##' @author Chris Wallace
##' @export
##' @keywords htest
##' @examples
##' 
##' data(for.exercise,package="snpStats")
##' case <- snps.10[subject.support$cc==1,]
##' control <- snps.10[subject.support$cc==0,]
##' summary(pairtest(case,control))
##' 
pairtest <- function(case,control,n.perm=0,pheno.perm=NULL,quiet=FALSE) {
  if(!is(case,"SnpMatrix") || !is(control,"SnpMatrix") || !identical(colnames(case),colnames(control))) {
    print(class(case))
    print(class(control))
    print(colnames(case))
    print(colnames(control))
    stop("case and control SnpMatrix objects must contain the same SNPs, in the same order\n")
  }
  d <- new("SnpMatrix",rbind2(case,control))
  if(!is.null(pheno.perm)) {
    if(!is.matrix(pheno.perm))
      pheno.perm <- as.matrix(pheno.perm, ncol=1)
    if(nrow(pheno.perm) != nrow(d))
      stop("length of phenotype vector must equal total number of cases and controls.\n")
    n.perm <- ncol(pheno.perm)
  }
  pheno <- rep(c(1,0),times=c(nrow(case),nrow(control)))
  permute <- n.perm > 0
  if(!permute) {
    return(p.value(single.snp.tests(phenotype=pheno,snp.data=d),df=1))
  }
  result <- matrix(NA,ncol(case),n.perm)
  if(is.null(pheno.perm))
    pheno.perm <- genperms(pheno,n.perm)
  for(i in 1:n.perm) {
    if(!quiet)
      cat(".")
    result[,i] <- p.value(single.snp.tests(phenotype=pheno.perm[,i],snp.data=d),df=1)
  }
  return(result)
}
