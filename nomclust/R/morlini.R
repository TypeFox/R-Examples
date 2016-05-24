#' Morlini and Zani's Measure S2
#' 
#' @description The S2 measure was proposed by Morlini and Zani (2012)
#' and it is based on a transformed dataset, which contains only binary variables (dummy coding).
#' Hierarchical clustering methods require a proximity (dissimilarity) matrix instead of a similarity matrix as
#' an entry for the analysis; therefore, dissimilarity \code{D} is computed from similarity \code{S} according the equation
#' \code{1/S-1}.\cr
#' \cr              
#' The use and evaluation of clustering with this measure can be found e.g. in (Sulc and Rezankova, 2014) or (Sulc, 2015).
#' \cr 
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
#' Morlini, I., Zani, S. (2012). A new class of weighted similarity indices using polytomous variables.
#' In Journal of Classification, 29(2), p. 199-226.
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
#' \code{\link[nomclust]{of}},
#' \code{\link[nomclust]{sm}}.
#'
#' @author Zdenek Sulc. \cr Contact: \email{zdenek.sulc@@vse.cz}
#' 
#' @examples
#' #sample data
#' data(data20)
#' # Creation of proximity matrix
#' prox_morlini <- morlini(data20)
#' 
#' @export 




morlini <- function(data) {
  
  #require(dummies)
  
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
  
  
  s <- ncol(data)
  num_cat <- sapply(data, function(x) length(unique(x)))
  
  #abs.freq <- freq.abs(data)
  
  data_dummy <- dummy.data.frame(data, dummy.classes ="ALL")
  
  n <- nrow(data_dummy)
  hs <- ncol(data_dummy)
  
  nsv <- sapply(data_dummy, sum)
  fsv2 <- log(1/(nsv/n)^2)
  
  E <- matrix(data=0,nrow=n,ncol=n)
  agreement <- vector(mode="numeric", length=hs)
  
  #computation of Eij
  for (i in 1:n) {
    for (j in 1:n) {
      for (k in 1:hs) {
        if (data_dummy[i,k] == 1 & data_dummy[j,k] == 1) {
          agreement[k] <- 1
        }
        else {
          agreement[k] <- 0
        }
      }
      E[i,j] <- fsv2 %*% agreement
    }
  }

#computation of Fij
  cum <- cumsum(num_cat)
  F <- matrix(data=0,nrow=n,ncol=n)

  for (i in 1:n) {
    for (j in 1:n) {
      v <- 0
      agreement <- vector(mode="numeric", length=hs)
      for (k in 1:s) {
        for (t in (v+1):cum[k]) {
          if (data_dummy[i,t] == 0 & data_dummy[j,t] == 1) {
            agreement[(v+1):cum[k]] <- 1
          }
        }
        v <- cum[k]      
      }
      F[i,j] <- t(agreement) %*% fsv2
    }
  }
  morlini <- 1 - E/(E+F)
  return(morlini)
}
  


  
  