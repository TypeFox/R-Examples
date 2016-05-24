#' Transform bases into genotypes
#' 
#' Preprocess the input data in case bases are provided instead of genotypes
#'  
#' @param input Matrix of genotype bases. Rows represent the individual markers. 
#' Columns represent samples, dependenden on the ploidy.
#' @param ploidy Ploidy level of the organism. 
#' Influences how many columns are collapsed into one.
#' @return Matrix of genotypes. 
#' The number of columns is 1/\code{ploidy} of the \code{input}.
#' @examples
#' data(simTetra)
#' bases2genotypes(simTetra, 4)
#' @export
bases2genotypes <- function(input,ploidy){
  nGeno<-ncol(input)
  output<-apply(input, 1, function(y) sapply(seq(1, nGeno, by = ploidy), 
                                             function(x) sum(y[x:(x + ploidy - 1)] == y[1])))
  t(output)
}

#' Randomize marker order and alleles within samples
#'  
#' In simulated datasets, the order or markers and alleles within samples is often given.
#' To remove any prior knowledge, that would not be available, the data should be randomized.
#' Thus, the performance of our tool can be validated unbiased.
#' 
#' @param input Matrix of genotypes. Rows represent markers. 
#' Columns represent samples.
#' @param ploidy Ploidy level of the organism. Default is 4.
#' @param ignore In case of unnecessary fronstanding columns (e.g. parental genotypes or rownames), these can be excluded from the randomization.
#' @return Matrix of the same size as the input matrix.
#' The markers are in a random order and the alleles within the samples are in a random order.
#' @examples
#' data(simTetra)
#' shuffleInput(simTetra, 4)
#' @export
shuffleInput <- function(input, ploidy = 4, ignore = 0){
  numbMark <- nrow(input)
  numbSamp <- ncol(input) / ploidy
  if(numbSamp %% 1 != 0){
    stop("Number of columns (", ncol(input) , ") is not ploidy (", ploidy,
         ") * number of samples (", numbSamp, ").")
  }
  #shuffle rows
  shufMat <- input[sample(1:numbMark, numbMark, replace = FALSE), ]
  #shuffle columns individual wise
  if(ploidy > 1){
    for(i in seq(ignore + 1, ncol(shufMat), by = ploidy)){
      for(j in 1:numbMark){
        vec <- i : (i + ploidy - 1)
        shufMat[j, vec] <- sample(shufMat[j, vec], ploidy, replace = FALSE)
      }
    }  
  }
  return(shufMat)
}

