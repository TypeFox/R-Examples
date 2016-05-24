#' get useful outputs from simulated annealing processes
#' @param simAnnealList the output from simAnnealing process
#' @param num number of simulated annealing runs
#' @return a list of three elements
#' 
#'  \item{costs_all}{costs of all simulated annealing runs.}
#'  \item{bestRankOrder}{best rank order found in all simulated annealing processes}
#'  \item{allRankOrder}{a dataframe, all best rank orders found in each simulated annealing processes}
getSimOutput <- function(simAnnealList, num){
  costs_all <- unlist(do.call(rbind, lapply(simAnnealList, function(x)x$Cb)))
  best_cost_index <- which.min(costs_all)
  allRankOrder <- lapply(simAnnealList, function(x) x$Ordb)
  names(allRankOrder) <- paste0("SimAnneal", 1:num)
  allRankOrder_df <- data.frame(allRankOrder)
  bestRankOrderIndex <- allRankOrder[[best_cost_index]]
  return(list(
    costs_all = costs_all,
    bestRankOrder = bestRankOrderIndex,
    allRankOrder = allRankOrder_df
  ))
}


#' Associate each costs with its corresponding simulated annealing runs
#' @param costs_all costs of all simAnnealing runs. It is the first element of the output from \code{getSimOutput}.
#' @param num number of simulated annealing runs
#' @return a data.frame of all costs.
getAllCosts <- function(costs_all, num){
  # export costs for each simAnnealRun
  CostOutput <- data.frame(simAnnealRun = c(1:num), Cost = costs_all)
  return(CostOutput)
}

#' assign IDs to the best rank order
#' @param ID_index it depends on the inputed data from simRankOrder. It takes the colnames of data as ID, and index this ID by its position in the colname.
#' @param bestRankOrder the best rank order found in all simulated annealing runs. It is the second output from \code{getSimOutput}.
#' @return a data.frame of all costs.
getBestRankOrder <- function(ID_index, bestRankOrder) {
  bestRankOrder <- data.frame(ID = ID_index[bestRankOrder, "ID"],
                              ranking = 1:nrow(ID_index))
  return(bestRankOrder)
}

#' assign IDs to all best rank orders
#' @param ID_index it depends on the inputed data from simRankOrder. It takes the colnames of data as ID, and index this ID by its position in the colname.
#' @param allRankOrder all rank orders found in all simulated annealing runs. It is the third output from \code{getSimOutput}.
#' @return a data.frame of all costs.
getAllRankOrder <- function(ID_index, allRankOrder){
  
  allRankOrder_df <- data.frame(apply(allRankOrder, 
                                      2, 
                                      function(x) ID_index[x, "ID"]),
                                stringsAsFactors = FALSE)
  return(allRankOrder_df)
}
