#' Goodall 4 Measure
#' 
#' @description The Goodall 4 similarity measure was firstly introduced in (Boriah et al., 2008).                                  
#' The measure ssigns higher similarity if the frequent categories match.
#' When measuring similarity between two variables, this measure provides complement results of Goodall 3 to one.
#' Hierarchical clustering methods require a proximity (dissimilarity) matrix instead of a similarity matrix as
#' an entry for the analysis; therefore, dissimilarity \code{D} is computed from similarity \code{S} according the equation
#' \code{1/S-1}.\cr
#' \cr                                    
#' The use and evaluation of clustering with this measure can be found e.g. in (Sulc, 2015).                                     
#'  
#' @param data data frame with cases in rows and variables in colums. Cases are characterized by nominal (categorical) variables coded as numbers.
#' 
#' @return Function returns a matrix of the size \code{n x n}, where \code{n} is the number of objects in original data. The matrix contains proximities
#' between all pairs of objects. It can be used in hierarchical cluster analyses (HCA), e.g. in \code{\link[cluster]{agnes}}.
#' \cr
#' @references
#' Boriah, S., Chandola and V., Kumar, V. (2008). Similarity measures for categorical data: A comparative evaluation.
#' In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254. Available at:
#' \url{ http://www-users.cs.umn.edu/~sboriah/PDFs/BoriahBCK2008.pdf}.
#' \cr
#' \cr
#' Goodall, V.D. (1966). A new similarity index based on probability. Biometrics. Vol. 22, No.4, p. 882.
#' \cr
#' \cr
#' Sulc, Z. (2015). Application of Goodall's and Lin's similarity measures in hierarchical clustering.
#' In Sbornik praci vedeckeho seminare doktorskeho studia FIS VSE. Praha: Oeconomica, 2015, p. 112-118. Available at:
#' \url{http://fis.vse.cz/wp-content/uploads/2015/01/DD_FIS_2015_CELY_SBORNIK.pdf}.
#'
#' @seealso
#' \code{\link[nomclust]{eskin}},
#' \code{\link[nomclust]{good1}},
#' \code{\link[nomclust]{good2}},
#' \code{\link[nomclust]{good3}},
#' \code{\link[nomclust]{iof}},
#' \code{\link[nomclust]{lin}},
#' \code{\link[nomclust]{lin1}},
#' \code{\link[nomclust]{morlini}}
#' \code{\link[nomclust]{of}},
#' \code{\link[nomclust]{sm}}.
#'
#' @author Zdenek Sulc. \cr Contact: \email{zdenek.sulc@@vse.cz}
#' 
#' @examples
#' #sample data
#' data(data20)
#' # Creation of proximity matrix
#' prox_goodall_4 <- good4(data20)
#'
#' @export

good4 <- function(data) {
  
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
  freq.rel <- freq.abs/r
  
  agreement <- vector(mode="numeric", length=s)
  good4 <- matrix(data=0,nrow=r,ncol=r)
  
  for (i in 1:r) {
    for (j in 1:r) {
      for (k in 1:s) {
        c <- data[i,k]
        if (data[i,k] == data[j,k]) {
          agreement[k] <- freq.rel[c,k]^2
        }
        else {
          agreement[k] <- 0
        }
      }
      if (i == j) {
        good4[i,j] <- 0
      }
      else {
        good4[i,j] <- 1-1/s*(sum(agreement))
      }
    }
  }
  return(good4)
}