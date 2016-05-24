#' Lin 1 Measure
#' 
#' @description The Lin 1 similarity measure was firstly introduced in (Boriah et al., 2008). In has
#' a complex system of weights. In case of mismatch, lower similarity is assigned if either
#' the mismatching values are very frequent or their relative frequency is in between the relative
#' frequencies of mismatching values. Higher similarity is assigned if the mismatched categories
#' are infrequent and there are a few other infrequent categories. In case of match,
#' lower similarity is given for matches on frequent categories or matches on categories
#' that have many other values of the same frequency. Higher similarity is given to matches
#' on infrequent categories.
#'                                        
#' Hierarchical clustering methods require a proximity (dissimilarity) matrix instead of a similarity matrix as
#' an entry for the analysis; therefore, dissimilarity \code{D} is computed from similarity \code{S} according the equation
#' \code{1/S-1}. After this transformation, it may
#' happen that some values in a proximity matrix get the value \code{-Inf}. Therefore, the following adjustment is applied:
#' \code{max(prox)+1}, where \code{prox} is a proximity matrix.
#' \cr
#' \cr
#'                                        
#' The use and evaluation of clustering with this measure can be found e.g. in (Sulc, 2015).
#'  
#' @param data data frame with cases in rows and variables in colums. Cases are characterized by nominal (categorical) variables coded as numbers.
#' 
#' @return Function returns a matrix of the size \code{n x n}, where \code{n} is the number of objects in original data. The matrix contains proximities
#' between all pairs of objects. It can be used in hierarchical cluster analyses (HCA), e.g. in \code{\link[cluster]{agnes}}.
#' \cr
#'
#' @references
#' Boriah, S., Chandola and V., Kumar, V. (2008). Similarity measures for categorical data: A comparative evaluation.
#'  In: Proceedings of the 8th SIAM International Conference on Data Mining, SIAM, p. 243-254. Available at:
#'  \url{ http://www-users.cs.umn.edu/~sboriah/PDFs/BoriahBCK2008.pdf}.
#'  \cr
#'  \cr
#' Sulc, Z. (2015). Application of Goodall's and Lin's similarity measures in hierarchical clustering.
#' In Sbornik praci vedeckeho seminare doktorskeho studia FIS VSE. Praha: Oeconomica, 2015, p. 112-118. Available at:
#' \url{http://fis.vse.cz/wp-content/uploads/2015/01/DD_FIS_2015_CELY_SBORNIK.pdf}.
#'
#' @seealso
#' \code{\link[nomclust]{eskin}},
#' \code{\link[nomclust]{good1}},
#' \code{\link[nomclust]{good2}},
#' \code{\link[nomclust]{good3}},
#' \code{\link[nomclust]{good4}},
#' \code{\link[nomclust]{iof}},
#' \code{\link[nomclust]{lin}},
#' \code{\link[nomclust]{morlini}},
#' \code{\link[nomclust]{of}},
#' \code{\link[nomclust]{sm}}.
#'
#' @author Zdenek Sulc. \cr Contact: \email{zdenek.sulc@@vse.cz}
#' 
#' @examples
#' #sample data
#' data(data20)
#' # Creation of proximity matrix
#' prox_lin1 <- lin1(data20)
#'
#' @export 


lin1 <- function(data) {
  
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
  freq.ln <- log(freq.rel)
  freq.ln[freq.ln == -Inf] <- 0
  
  
  agreement <- vector(mode="numeric", length=s)
  lin1 <- matrix(data=0,nrow=r,ncol=r)
  weights <- vector(mode="numeric", length=s)

  for (i in 1:r) {
    for (j in 1:r) {
      for (k in 1:s) {
        c <- data[i,k]
        d <- data[j,k]
        if (data[i,k] == data[j,k]) {
          logic <- freq.rel[,k] == freq.rel[c,k]
          agreement[k] <- sum(logic * freq.ln[,k])
          weights[k] <- sum(logic * freq.ln[,k])
        }
        else {
          if (freq.rel[c,k] >= freq.rel[d,k]) {
            logic <- freq.rel[,k] >= freq.rel[d,k] & freq.rel[,k] <= freq.rel[c,k]
            agreement[k] <- 2*log(sum(logic * freq.rel[,k]))
            weights[k] <- sum(logic * freq.ln[,k])
          }
          else {
            logic <- freq.rel[,k] >= freq.rel[c,k] & freq.rel[,k] <= freq.rel[d,k]
            agreement[k] <- 2*log(sum(logic * freq.rel[,k]))
            weights[k] <- sum(logic * freq.ln[,k])
          }
        }
      }
      lin1[i,j] <- 1/(1/sum(weights)*(sum(agreement))) - 1
    }
  }
  lin1[lin1 == -Inf] <- max(lin1) + 1
  return(lin1)
}