#' Implementation of MULTIMOORA Method for Multi-Criteria Decision Making Problems.
#'
#' @description The \code{MMOORA} function implements both the Multi-Objetive Optimization by Ration Analysis (MOORA) and the "Full Multiplicative Form" (MULTIMOORA).
#' @param decision The decision matrix (\emph{m} x \emph{n}) with the values of the \emph{m} alternatives, for the \emph{n} criteria. 
#' @param weights A vector of length \emph{n}, containing the weights for the criteria. The sum of the weights has to be 1.
#' @param cb A vector of length \emph{n}. Each component is either \code{cb(i)='max'} if the \emph{i-th} criterion is benefit or \code{cb(i)='min'} if the \emph{i-th} criterion is a cost.
#' @return \code{MMOORA} returns a data frame which contains the scores and the four rankings calculated (Ratio System, Reference Point, Multiplicative Form and Multi-MOORA ranking). 
#' @references Brauers, W. K. M.; Zavadskas, E. K. Project management by MULTIMOORA as an instrument for transition economies. Technological and Economic Development of Economy, 16(1), 5-24, 2010.
#' @examples
#' 
#'  d <- matrix(rpois(12, 5), nrow = 4)
#'  w <- c(0.2, 0.2, 0.6)
#'  cb <- c('max','min','max')
#'  MMOORA(d,w,cb)

MMOORA <- function(decision, #matrix with all the alternatives
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
  
  
  #MMOORA method
  
  #1. Normalization and weighting
  d = sqrt(colSums(decision^2))
  NW <- matrix(nrow = nrow(decision), ncol = ncol(decision))
  for(j in 1:ncol(decision)){
    NW[,j] <- (decision[,j] / d[j]) * weights[j]
  }

  #2. Ration system
  NR <- NW
  for(j in 1:ncol(decision)){
    if (cb[j] == 'min'){
      NR[,j] <- NW[,j]*(-1)    
    }
  }
  RS <- apply(NR, 1, sum)
  
  #3. Reference point
  Ref <- as.integer(cb == "max") * apply(NW, 2, max) + 
    as.integer(cb == "min") * apply(NW, 2, min)
  RefP <- matrix(nrow = nrow(decision), ncol = ncol(decision))
  for(j in 1:ncol(decision)){
    RefP[,j] <- abs(Ref[j]-NW[,j])
  }
  RP <- apply(RefP, 1, max)
  
  #4. Multiplicative form
  NEW <- matrix(nrow = nrow(decision), ncol = ncol(decision))
  for(j in 1:ncol(decision)){
    NEW[,j] <- (decision[,j] / d[j]) ^ weights[j]
  }
  
  max <- NEW
  min <- NEW
  for (j in 1:ncol(NEW)){
    if (cb[j] == 'max'){
      min[,j] <- 1
    }else{
      max[,j] <- 1
    } 
  }
  A <- apply(max, 1, prod)
  B <- apply(min, 1, prod)
  M <- A/B

  #5. Ranking the alternatives
  Rrs <- rank(-RS)
  Rrp <- rank(RP)
  Rm <- rank(-M)
  R <- Rrs + Rrp + Rm
  return(data.frame(Alternatives = 1:nrow(decision), RatioSystem = RS, Ranking = rank(-RS, ties.method= "random"), ReferencePoint = RP, Ranking = rank(RP, ties.method= "random"), MultiplicativeForm = M, Ranking = rank(-M, ties.method= "random"), MultiMooraRanking = rank(R, ties.method= "random")))
  
}