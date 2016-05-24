#' Identify outliers in a variable
#' 
#'  Takes in a vector, and returns count and index of outliers
#'  @param vector an integer or numeric vector
#'  @details
#'  The function uses the same criteria to identify outliers as the one used for
#'  box plots. All values that are greater than 75th percentile value + 1.5 times the 
#'  inter quartile range or lesser than 25th percentile value - 1.5 times the inter
#'  quartile range, are tagged as outliers.
#'  
#'  The individual elements (number of outliers and index of outliers) of the two 
#'  element output list can be picked using the code given in example. The index 
#'  of outliers can be used to get a vector of all outliers.
#'  @return a list with two elements: count and index of outliers
#'  @author Akash Jain
#'  @seealso \code{\link{decile}}, \code{\link{pentile}}, \code{\link{imputemiss}}
#'  @examples
#'  # Scores vector
#' scores <- c(1, 4, 7, 10, 566, 21, 25, 27, 32, 35, 
#'             49, 60, 75, 23, 45, 86, 26, 38, 34, 223, -3)
#' 
#' # Identify the count of outliers and their index
#' ltOutliers <- outliers(vector = scores)
#' numOutliers <- ltOutliers$numOutliers
#' idxOutliers <- ltOutliers$idxOutliers
#' valOutliers <- scores[idxOutliers]
#'  @export
outliers <- function(vector) {
  if(class(vector) != 'integer' && class(vector) != 'numeric') {
    stop('Invalid input: vector should be either integer or numeric')
  } else {
    p25 <- quantile(vector, c(0.25), na.rm = TRUE)
    p75 <- quantile(vector, c(0.75), na.rm = TRUE)
    iqr <- p75 - p25
    uplim <- p75 + 1.5*iqr
    lowlim <- p25 - 1.5*iqr
    numOutliers <- sum(vector < lowlim, na.rm = TRUE) + sum(vector > uplim, na.rm = TRUE)
    idxOutliers <- sort(c(match(na.omit(vector[vector < lowlim]), vector), 
                          match(na.omit(vector[vector > uplim]), vector)))
    valueOutliers <- vector[idxOutliers]
    return(list(numOutliers = numOutliers, 
                idxOutliers = idxOutliers))    
  }
}