#' Calculate the Damerau Levenshtein Distance between two strings
#'
#' The restricted Damerau Levenshtein Distance between two strings is commonly used for checking typographical errors in strings.
#' It takes the deletion and insertion of a character, a wrong character (substition) or the swapping (transposition) 
#' of two characters into account. By default these operations each account for distance 1.
#'
#' @references
#' 
#' Damerau F (1964). A technique for computer detection and correction of
#' spelling errors. Communications of the ACM, 7,issue 3
#'
#' Levenshtein VI (1966). Binary codes capable of correcting deletions, insertions, 
#' and reversals. Soviet Physics Doklady 10: 707-10
#' Damerau Levenshtein Distance calculates the difference between two strings
#' used for typographical errors (typo's) 
#'
#' @param sa character vector
#' @param sb character vector of equal \code{length(sa)}
#' @param w integer vector for cost of deletion, insertion, substitution 
#' and transposition.
#' 
#' @return integer vector with pairwise edit distances
damerauLevenshteinDistance <- function(sa,sb, w=c(1,1,1,1)){
   if (length(sa) != length(sb)) stop('sa and sb must be of equal length')
   mapply(function(a,b){      
      a <- c("",a)
      b <- c("",b)
      
      la <- length(a)
      lb <- length(b)
      d <- matrix(0, nrow = la, ncol = lb)
    
      
      d[,1] <- 0:(la-1)
      d[1,] <- 0:(lb-1)
      
      eq <- outer(a,b, "==")
      
      for(i in 2:la){
         for(j in 2:lb) {
         
            if (eq[i,j]){
               d[i,j] <- d[i-1,j-1]
               next
            }

            cost <- c( d[i-1, j  ]      # deletion
                     , d[i  , j-1]      # insert
                     , d[i-1, j-1]      # substitution
                     ) + w[1:3]

            if (i>2 & j>2 & eq[i-1, j] & eq[i, j-1])
                cost <- c(cost, d[i-2,j-2] + w[4]) 

            d[i,j] <- min(cost) 
         }
      }
      d[la,lb]
   }
   , strsplit(as.character(sa), NULL)
   , strsplit(as.character(sb), NULL)
   )
}
