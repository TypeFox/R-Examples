#' Occurence Frequency (OF) Measure
#' 
#' @description The OF (Occurrence Frequency) measure was originally constructed for the text mining,
#' see (Sparck-Jones, 1972), later, it was adjusted for categorical variables.
#' It assigns higher similarity to mismatches on less frequent values and otherwise. 
#' Hierarchical clustering methods require a proximity (dissimilarity) matrix instead of a similarity matrix as
#' an entry for the analysis; therefore, dissimilarity \code{D} is computed from similarity \code{S} according the equation
#' \code{1/S-1}.\cr
#' \cr                          
#'                                        
#'  
#' @param data data frame or matrix with cases in rows and variables in colums. Cases are characterized by nominal (categorical) variables coded as numbers.
#' 
#' @return Function returns a matrix of the size \code{n x n}, where \code{n} is the number of objects in original data. The matrix contains proximities
#' between all pairs of objects. It can be used in hierarchical cluster analyses (HCA), e.g. in \code{\link[cluster]{agnes}}.
#' \cr
#'
#' @references
#' Boriah, S., Chandola and V., Kumar, V. (2008). Similarity measures for categorical data: A comparative evaluation.
#' In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254. Available at:
#'  \url{ http://www-users.cs.umn.edu/~sboriah/PDFs/BoriahBCK2008.pdf}.
#'  \cr
#'  \cr
#' Spark-Jones, K. (1972). A statistical interpretation of term specificity and its application in retrieval.
#' In Journal of Documentation, 28(1), p. 11-21. Later: Journal of Documentation, 60(5) (2002), p. 493-502.
#' \cr
#' \cr
#' Sulc, Z. and Rezankova, H. (2014). Evaluation of recent similarity measures for categorical data.
#' In: AMSE. Wroclaw: Wydawnictwo Uniwersytetu Ekonomicznego we Wroclawiu, p. 249-258.
#' Available at: \url{http://www.amse.ue.wroc.pl/papers/Sulc,Rezankova.pdf}.
#' 
#' @seealso
#' \code{\link[nomclust]{eskin}},
#' \code{\link[nomclust]{good1}},
#' \code{\link[nomclust]{good2}},
#' \code{\link[nomclust]{good3}},
#' \code{\link[nomclust]{good4}},
#' \code{\link[nomclust]{iof}},
#' \code{\link[nomclust]{lin}},
#' \code{\link[nomclust]{lin1}},
#' \code{\link[nomclust]{morlini}},
#' \code{\link[nomclust]{sm}}.
#'
#' @author Zdenek Sulc. \cr Contact: \email{zdenek.sulc@@vse.cz}
#' 
#' @examples
#' #sample data
#' data(data20)
#' # Creation of proximity matrix
#' prox_of <- of(data20)
#' 
#' @export 




of <- function(data) {
  
  r <- nrow(data)
  s <- ncol(data)
  
  #recoding variables
  num_var <- ncol(data)
  num_row <- nrow(data)
  data2 <- matrix(data = 0, nrow = num_row, ncol = num_var)
  for (k in 1:num_var) {
    categories <- unique(data[, k])
    cat_new <- 1:length(categories)
    for (l in 1:length(categories)) {
      for (i in 1:num_row) {
        if (data[i, k] == categories[l]) {
          data2[i, k] <- cat_new[l]
        }
      }
    }
  }
  data <- data.frame(data2)
  
  
  
  freq.abs <- freq.abs(data)
  
  agreement <- vector(mode="numeric", length=s)
  of <- matrix(data=0,nrow=r,ncol=r)
  
  for (i in 1:r) {
    for (j in 1:r) {
      for (k in 1:s) {
        c <- data[i,k]
        d <- data[j,k]
        if (data[i,k] == data[j,k]) {
          agreement[k] <- 1
        }
        else {
          agreement[k] <- 1/(1+log(r/freq.abs[c,k])*log(r/freq.abs[d,k]))
        }
      }
      of[i,j] <- 1/(1/s*(sum(agreement)))-1
    }
  }
  return(of)
}

