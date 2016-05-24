#' Random generation of odds ratio (OR) effect sizes.  
#'
#' Generates random odds ratios, logged odds ratios, and their 
#' variances (Cornfield 1951).       
#'
#' @param K Number of effect sizes to generate.
#' @param p_A The odds of the event of interest for Group A.  A probability 
#'    ranging from zero to one.  
#' @param p_B The odds of the event of interest for Group B.  A probability 
#'    ranging from zero to one. 
#' @param N_A The total number of samples of Group A.
#' @param N_B The total number of samples of Group B.
#' @param continuity Odds ratios with zero events cannot be computed.  Following, 
#'    Cox (1970), a continuity correction can be added to each cell of the 2 by 2
#'    table to help improve this problem of zero events within the table.  The
#'    default value added is 0.5.
#' @param logged When \code{"FALSE"}, returns non-logged transformed 
#'    odds ratios and appropriate variances.  Default is TRUE.
#'
#' @return A data table with columns of random effect sizes (OR) and their 
#'    variances.
#'
#' @examples
#'    random_OR(K = 5, p_A = 0.3, N_A = 100, p_B = 0.1, N_B = 60)
#'
#' @references Cornfield, J. 1951. A method for estimating comparative rates
#'    from Clinical Data. Applications to cancer of the lung, breast, and cervix.
#'    Journal of the National Cancer Institute 11: 1269-1275.
#' @references Cox, D.R. 1970. The continuity correction. Biometrika 57: 217-219.
#'
#' @importFrom stats rmultinom
#' @export random_OR

random_OR <- function(K,
                      p_A,
                      N_A,
                      p_B,
                      N_B,
                      continuity = 0.5,
                      logged = TRUE) {
   
  outcomes_A <- rmultinom(K , size = N_A, prob = c(p_A, 1.0 - p_A))
  outcomes_B <- rmultinom(K , size = N_B, prob = c(p_B, 1.0 - p_B))

  continuityFix <- apply(rbind(outcomes_A, outcomes_B), 2, function(x) 0 %in% x)
  outcomes_A[, continuityFix] <- outcomes_A[, continuityFix] + continuity
  outcomes_B[, continuityFix] <- outcomes_B[, continuityFix] + continuity
  
  theOdds <- (outcomes_A[1,]/outcomes_B[1,])/(outcomes_A[2,]/outcomes_B[2,])
  var_theOdds <- 1.0/outcomes_A[1,] + 1.0/outcomes_B[1,] + 1.0/outcomes_A[2,] + 1.0/outcomes_B[2,]
  
  if(logged == TRUE) {
    theOdds <- log(theOdds)
  } else {
    var_theOdds <- (theOdds ^ 2.0) * var_theOdds
  }
  
  random_OR <- data.frame(OR = theOdds, var_OR = var_theOdds)
  return(random_OR)
}