#' Implementation of TOPSIS Method for Multi-Criteria Decision Making Problems.
#'
#' @description The \code{TOPSIS} function implements the Technique for Order of Preference by Similarity to Ideal Solution (TOPSIS) Method.
#' @param decision The decision matrix (\emph{m} x \emph{n}) with the values of the \emph{m} alternatives, for the \emph{n} criteria.  
#' @param weights A vector of length \emph{n}, containing the weights for the criteria. The sum of the weights has to be 1.
#' @param cb A vector of length \emph{n}. Each component is either \code{cb(i)='max'} if the \emph{i-th} criterion is benefit or \code{cb(i)='min'} if the \emph{i-th} criterion is a cost.
#' @return \code{TOPSIS} returns a data frame which contains the score of the R index and the ranking of the alternatives.
#' @references Hwang, C.L.; Yoon, K. Multiple Attribute Decision Making. In: Lecture Notes in Economics and Mathematical Systems 186. Springer-Verlag, Berlin, 1981.
#' @examples
#' 
#'  d <- matrix(rpois(12, 5), nrow = 4)
#'  w <- c(0.2, 0.2, 0.6)
#'  cb <- c('max','min','max')
#'  TOPSIS(d,w,cb)

TOPSIS <- function(decision, #matrix with all the alternatives
                   weights,  #vector with the numeric values of the weights
                   cb        #vector with the "type" of the criteria (benefit = "max", cost = "min")
)
{
  #Checking the arguments
  if(! is.matrix(decision))
    stop("'decision' must be a matrix with the values of the alternatives")
  if(missing(weights))
    stop("a vector containing n weigths, adding up to 1, should be provided")
  if(sum(weights) != 1)
    stop("The sum of 'weights' is not equal to 1")
  if(! is.character(cb))
    stop("'cb' must be a character vector with the type of the criteria")
  if(! all(cb == "max" | cb == "min"))
    stop("'cb' should contain only 'max' or 'min'")
  if(length(weights) != ncol(decision))
    stop("length of 'weights' does not match the number of the criteria")
  if(length(cb) != ncol(decision))
    stop("length of 'cb' does not match the number of the criteria")
  
  
  #TOPSIS method
  
  #1. Normalization and weighting
  d = sqrt(colSums(decision^2))
  NW <- matrix(nrow = nrow(decision), ncol = ncol(decision))
  for(j in 1:ncol(decision)){
    NW[,j] <- (decision[,j] / d[j]) * weights[j]
  }

  #2. Ideal solutions
  posI <- as.integer(cb == "max") * apply(NW, 2, max) + 
    as.integer(cb == "min") * apply(NW, 2, min)
  negI <- as.integer(cb == "min") * apply(NW, 2, max) + 
    as.integer(cb == "max") * apply(NW, 2, min)
  
  #3. Distances to the ideal solutions
  distance =function(x,y){
    sqrt(sum((x - y) ^ 2))
  }
  posDis <- apply(NW, 1, distance, posI)
  negDis <- apply(NW, 1, distance, negI)
  
  #4. R index
  R <- negDis/(negDis+posDis)
  
  #5. Rank the alternatives
  return(data.frame(Alternatives = 1:nrow(decision), R = R, Ranking = rank(-R, ties.method= "random")))
  
}