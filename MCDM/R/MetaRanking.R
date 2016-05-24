#' Implementation of MetaRanking function for Multi-Criteria Decision Making Problems.
#'
#' @description The \code{MetaRanking} function internally calls functions  \code{MMOORA}, \code{TOPSIS}, \code{VIKOR} and \code{WASPAS} and then calculates a sum of the their rankings.
#' @param decision The decision matrix (\emph{m} x \emph{n}) with the values of the \emph{m} alternatives, for the \emph{n} criteria. 
#' @param weights A vector of length \emph{n}, containing the weights for the criteria. The sum of the weights has to be 1.  
#' @param cb A vector of length \emph{n}. Each component is either \code{cb(i)='max'} if the \emph{i-th} criterion is benefit or \code{cb(i)='min'} if the \emph{i-th} criterion is a cost.
#' @param lambda A value in [0,1]. It is used in the calculation of the W index for WASPAS method.
#' @param v A value in [0,1]. It is used in the calculation of the Q index for VIKOR method.
#' @return \code{MetaRanking} returns a data frame which contains the rankings of the Multi-MOORA, TOPSIS, VIKOR, WASPAS Methods and the MetaRanking of the alternatives. 
#' @examples
#' 
#'  d <- matrix(rpois(12, 5), nrow = 4)
#'  w <- c(0.2, 0.2, 0.6)
#'  cb <- c('max','min','max')
#'  lambda <- 0.5
#'  v <- 0.5
#'  MetaRanking(d,w,cb,lambda,v)

MetaRanking <- function(decision, #matrix with all the alternatives
                        weights,  #vector with the numeric values of the weights
                        cb,       #vector with the "type" of the criteria (benefit = "max", cost = "min")
                        lambda,   #value with the real number of the 'lambda' parameter to calculate W
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
  if(missing(lambda))
    stop("a value for 'lambda' in [0,1] should be provided")
  if(missing(v))
    stop("a value for 'v' in [0,1] should be provided")
  
  
  #Multi-MOORA method 
  MMoora = MMOORA(decision,weights,cb)
  
  #TOPSIS method 
  Topsis = TOPSIS(decision,weights,cb)
  
  #VIKOR method  
  Vikor = VIKOR(decision,weights,cb,v)
  
  #WASPAS method  
  Waspas = WASPAS(decision,weights,cb,lambda)
  
  #Meta-Ranking
  MetaR = MMoora[,8]+Topsis[,3]+Vikor[,5]+Waspas[,3]
  return(data.frame(Alternatives = 1:nrow(decision), MMOORA = MMoora[,8], TOPSIS = Topsis[,3], VIKOR = Vikor[,5], WASPAS = Waspas[,3], METARANKING = rank(MetaR, ties.method= "random")))
  
}
