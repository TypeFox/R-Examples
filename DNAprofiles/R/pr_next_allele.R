#' Probability of seeing next allele (Dirichlet sampling)
#'
#' @param i Integer (vector), allele number
#' @param seen Integer matrix with alleles already seen
#' @param fr Numeric vector with allelic proportions
#' @param theta Numeric giving the inbreeding coefficient
#' @details When a population is subdivided into subpopulations, consecutively sampled alleles are not independent draws. This function implements the Dirichlet formula which states that after sampling \eqn{n} alleles, of which \eqn{m} are of type \eqn{A_i}, the probability that the next allele is of type \eqn{A_i} equals:
#'
#' \eqn{(m*\theta+(1-\theta)*p_i)/(1+(n-1)*\theta)}
#' 
#' The alleles already sampled are passed as the rows of the matrix \code{seen}, while the corresponding element of \code{i} specifies for which next allele the probability of sampling is computed. The length of \code{i} has to be equal to the number of rows of \code{seen}.
#' 
#' @return numeric (vector) of probabilities
#' @seealso \code{\link{pr.next.alleles}}, \code{\link{rmp}}
#' @examples
#' # theta=0 means independent sampling, so after seeing
#' # allele 1 three times the pr. remains 1/2
#' pr.next.allele(1,seen=matrix(c(1,1,1),nrow=1),fr=c(1/2,1/2),theta=0)
#' # theta>0 slighly increases the pr. of
#' # seeing the same allele again
#' pr.next.allele(1,seen=matrix(c(1,1,1),nrow=1),fr=c(1/2,1/2),theta=0.05)
#'
#' # the function also works on vectors
#' # after seeing 1,1,1, the pr. of 1 remains 1/2 
#' # and the same applies to the pr. of 2 after 2,2,1
#' pr.next.allele(c(1,2),seen=matrix(c(1,1,1,2,2,1),nrow=2,byrow=TRUE),fr=c(1/2,1/2),theta=0)
#' pr.next.allele(c(1,2),seen=matrix(c(1,1,1,2,2,1),nrow=2,byrow=TRUE),fr=c(1/2,1/2),theta=0.05)
#' @export
pr.next.allele <- function(i,seen,fr,theta=0){
  if (!is.matrix(seen)) stop("seen must be a matrix with n (the number of alleles) columns")
  if (length(i)!=nrow(seen)) stop("The length of i should be equal to the number of rows of seen")
  if (min(i,na.rm = TRUE)<1L) stop("ij should contain positive integers")
  if (max(i,na.rm = TRUE)>length(fr)) stop("ij cannot contain alleles outside of allele ladder fr") 
  if (any(is.na(seen))) stop("seen should not contain NAs")  
  
  n <- ncol(seen) # n total alleles seen
  m <- rowSums(matrix(apply(seen,2,function(s0) s0==i),nrow=nrow(seen)))# of which m are of type i    
  (m*theta+(1-theta)*fr[i])/(1+(n-1)*theta) # dirichlet formula
}
NULL
#' Probability of seeing next alleles (Dirichlet sampling)
#'
#' @param ij integer matrix with allele numbers
#' @param seen integer matrix with alleles already seen
#' @param fr numeric vector with allelic proportions
#' @param theta numeric background relatedness
#' @details When a population is subdivided into subpopulations, consecutively sampled alleles are not independent draws. This function implements the Dirichlet formula which states that after sampling \eqn{n} alleles, of which \eqn{m} are of type \eqn{A_i}, the probability that the next allele is of type \eqn{A_i} equals:
#'
#' \eqn{(m*\theta+(1-\theta)*p_i)/(1+(n-1)*\theta)}
#' 
#' The alleles already sampled are passed as the rows of the matrix \code{seen}, while the corresponding row of \code{ij} specifies for which alleles the probability of sampling is computed. The numer of rows of \code{ij} has to be equal to the number of rows of \code{seen}.
#' @return numeric (vector) of probabilities
#' @seealso \code{\link{pr.next.allele}}, \code{\link{rmp}}
#' @examples
#' # compute the pr. of seeing 1,1 after 1,1,1,1
#' # when theta=0 this is simply p_1^2
#' pr.next.alleles(t(c(1,1)),seen=t(c(1,1,1,1)),fr=c(1/4,3/4),theta=0)
#' # but when theta>0, the pr. of seeing more 1's increases slightly
#' pr.next.alleles(t(c(1,1)),seen=t(c(1,1,1,1)),fr=c(1/4,3/4),theta=0.05)
#' 
#' # pr. distribution of (ordered!) genotypes after seeing 1,1,1
#' ij=matrix(c(1,1,1,2,2,2),ncol=2,byrow=TRUE)
#' seen=matrix(1,nrow=3,ncol=3,byrow=TRUE)
#' pr.next.alleles(ij,seen,fr=c(1/4,3/4),theta=0) # theta=0
#' pr.next.alleles(ij,seen,fr=c(1/4,3/4),theta=0.1) # theta=0.1
#' 
#' p0 <- pr.next.alleles(ij,seen,fr=c(1/4,3/4),theta=0)
#' stopifnot(all.equal(p0[1]+2*p0[2]+p0[3],1))
#' 
#' p1 <- pr.next.alleles(ij,seen,fr=c(1/4,3/4),theta=0.05)
#' stopifnot(all.equal(p1[1]+2*p1[2]+p1[3],1))
#' @export
pr.next.alleles <- function(ij,seen,fr,theta=0){
  if (!is.matrix(seen)) stop("seen must be a matrix with n (the number of alleles already sampled) columns")
  if (!is.matrix(ij)) stop("ij must be matrix with at least 1 column")  
  if ((nrow(ij))!=nrow(seen)) stop("ij must contain as many rows as seen")
  if (min(ij,na.rm = TRUE)<1L) stop("ij should contain positive integers")
  if (max(ij,na.rm = TRUE)>length(fr)) stop("ij cannot contain alleles outside of allele ladder fr") 
  if (any(is.na(seen))) stop("seen should not contain NAs")
  if (any(is.na(ij))) stop("ij should not contain NAs")
  
  if (ncol(ij)>1){
    pr.next.allele(i=ij[,ncol(ij)],seen=seen,fr=fr,theta=theta)*
      pr.next.alleles(ij=ij[,-ncol(ij),drop=FALSE],seen=cbind(ij[,ncol(ij)],seen),fr=fr,theta=theta)
  }else{
    pr.next.allele(i=ij,seen=seen,fr=fr,theta=theta)
  }
}