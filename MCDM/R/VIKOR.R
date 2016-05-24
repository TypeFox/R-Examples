#' Implementation of VIKOR Method for Multi-Criteria Decision Making Problems.
#'
#' @description The \code{VIKOR} function implements the "VIseKriterijumska Optimizacija I Kompromisno Resenje" (VIKOR) Method.
#' @param decision The decision matrix (\emph{m} x \emph{n}) with the values of the \emph{m} alternatives, for the \emph{n} criteria. 
#' @param weights A vector of length \emph{n}, containing the weights for the criteria. The sum of the weights has to be 1.
#' @param cb A vector of length \emph{n}. Each component is either \code{cb(i)='max'} if the \emph{i-th} criterion is benefit or \code{cb(i)='min'} if the \emph{i-th} criterion is a cost.
#' @param v A value in [0,1]. It is used in the calculation of the Q index.
#' @return \code{VIKOR} returns a data frame which contains the score of the S, R and Q indixes and the ranking of the alternatives according to Q index. 
#' @references Opricovic, S.; Tzeng, G.H. Compromise solution by MCDM methods: A comparative analysis of VIKOR and TOPSIS. European Journal of Operational Research, 156(2), 445-455, 2004.
#' @examples
#' 
#'  d <- matrix(rpois(12, 5), nrow = 4)
#'  w <- c(0.2, 0.2, 0.6)
#'  cb <- c('max','min','max')
#'  v <- 0.5
#'  VIKOR(d,w,cb,v)

VIKOR <- function(decision, #matrix with all the alternatives
                  weights,  #vector with the numeric values of the weights
                  cb,       #vector with the "type" of the criteria (benefit = "max", cost = "min")
                  v         #value with the real number of the 'v' parameter to calculate Q
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
  if(missing(v))
    stop("a value for 'v' in [0,1] should be provided")
  
  
  #VIKOR method
  
  #1. Ideal solutions
  posI <- as.integer(cb == "max") * apply(decision, 2, max) + 
    as.integer(cb == "min") * apply(decision, 2, min)
  negI <- as.integer(cb == "min") * apply(decision, 2, max) + 
    as.integer(cb == "max") * apply(decision, 2, min)
  
  #2. S and R index
  norm =function(x,w,p,n){
    w*((p-x)/(p-n))
  }
  SAux <- apply(decision, 1, norm, weights, posI, negI)
  S <- apply(SAux, 2, sum)
  R <- apply(SAux, 2, max)
  
  
  #3. Q index
  Q <- v*(S-min(S))/(max(S)-min(S))+(1-v)*(R-min(R))/(max(R)-min(R))
  
  #4. Ranking the alternatives
  return(data.frame(Alternatives = 1:nrow(decision), S = S, R = R, Q = Q, Ranking = rank(Q, ties.method= "random")))
  
}