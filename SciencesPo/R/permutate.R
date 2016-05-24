#' @encoding UTF-8
#' @title  Create k random permutations of a vector
#' @description Creates a k random permutation of a vector.
#' @details should be used only for length(input)! >> k
#' @param input A vector to be permutated.
#' @param k number of permutations to be conducted.
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}.
#' @keywords Sampling, Manipulation
#' @examples
#' # row wise permutations
#' permutate(input=1:5, k=5)
#' @export
permutate <- function(input,k){
  n <- length(input)
  mat <- matrix(data=NA,nrow=k,ncol=n) # allocate memory
  k <- min(k, .npermutate(input))
  inserted <- 0
  while(inserted < k){
    p <- sample(input)
    # check if the vector has already been inserted
    if(sum(apply(mat,1,identical,p)) == 0){
      mat[inserted+1,] <- p
      inserted <- inserted+1
    }
  }
  return(mat)
}
NULL

.npermutate <- function(vec){
  tab <- table(vec); # count occurences of each element
  occurences <- tab[tab>1]; # get those greater than 1
  numerator <- lfactorial(length(vec))
  if(length(occurences ) > 0){
    denominator <- sum(sapply(occurences , lfactorial))
  } else {
    denominator <- 0
  }
 return(exp(numerator-denominator))
}
NULL
