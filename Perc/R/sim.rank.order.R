#' Find rank order using simulated annealing
#' 
#' \code{simRankOrder} find the rank order for the win-loss relationship
#' 
#' @param data a matrix. the win-loss probability matrix 
#' which is the second element of the output from \code{conductance}
#' @param num number of SimAnnealing (default is set at 10)
#' @param kmax an integer between 2 to 1000, indicating the number of simulations in each SimAnnealing.
#' @param alpha a positive integer that 
#' reflects the influence of an observed win/loss interaction 
#' on an underlying win-loss probability. 
#' It is used in the calculation of the posterior distribution 
#' for the win-loss probability of \code{i} over \code{j}: \eqn{Beta(\alpha c_{i,j} +\beta, c_{i,j}+\beta)}{Beta*(\alpha * c_ij + \beta, c_ij + \beta)}. 
#' In the absence of expertise to accurately estimate alpha, 
#' it is estimated from the data.
#' @return a list of two dataframes. 
#'    \item{BestSimulatedRankOrder}{a dataframe representing the best simulated rank order.}
#'    \item{Costs}{the cost of each simulated annealing run}
#'    \item{AllSimulatedRankOrder}{a dataframe representing all simulated rank orders.}
#' 
#' 
#' @references Fushing, H., McAssey, M. P., Beisner, B., & McCowan, B. (2011). Ranking network of a captive rhesus macaque society: a sophisticated corporative kingdom. PLoS One, 6(3), e17817-e17817.
#' 
#' @seealso \code{\link{conductance}} \code{\link{transitivity}}
#' 
#' @examples
#' # convert an edgelist to conflict matrix
#' confmatrix <- as.conflictmat(sampleEdgelist)
#' # find dominance probability matrix
#' perm2 <- conductance(confmatrix, maxLength = 2)
#' \dontrun{
#' # Note: It takes a while to run the simRankOrder example.
#' s.rank <- simRankOrder(perm2$p.hat, num = 10, kmax = 1000)
#' s.rank$BestSimulatedRankOrder
#' s.rank$Costs
#' s.rank$AllSimulatedRankOrder
#' }
#' \dontshow{
#' s.rank <- simRankOrder(perm2$p.hat, num = 2, kmax = 5)
#' s.rank$BestSimulatedRankOrder
#' s.rank$Costs
#' s.rank$AllSimulatedRankOrder
#' }
#' @export

simRankOrder <- function(data, num = 10, alpha = NULL, kmax = 1000){  # if null, take transitivity; if not null take specify
  # the output from percolation function.
  if (!(is.matrix(data))) {
    stop("The second element 'p.hat' from the output of 'conductance' should be used.")
  }
  
  if (any(data < 0)){
    stop("Values smaller than 0 detected. Please check your data. The second element 'p.hat' from the output of 'conductance' should be used.")
  }
  
  if (any(data > 1)){
    stop("Values greater than 1 detected. Please check your data. The second element 'p.hat' from the output of 'conductance' should be used.")
  }
  

  ### Run the simulated annealing many times, since it sometimes gets 
  ### stuck in local minima.
  
  # automatically replicate simAnneal for times user specified.
  simAnnealList <- replicate(num, SimAnneal(data, kmax), simplify = FALSE)
  
  ID_index <- data.frame(ID = colnames(data), 
                         index = 1:ncol(data),
                         stringsAsFactors = FALSE)
  simOutputs <- getSimOutput(simAnnealList = simAnnealList, num = num)
  bestRankOrder <- getBestRankOrder(ID_index = ID_index,
                                   bestRankOrder = simOutputs$bestRankOrder)
  CostOutput <- getAllCosts(simOutputs$costs_all, num = num)
  allRankingDF <- getAllRankOrder(ID_index = ID_index, 
                                  allRankOrder = simOutputs$allRankOrder)
  
  return(
    list(BestSimulatedRankOrder = bestRankOrder, 
         Costs = CostOutput, 
         AllSimulatedRankOrder = allRankingDF)
    )
}