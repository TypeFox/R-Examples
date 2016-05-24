#' Formats the estimated abundances of all species categories, to be consistent
#' with the prorated results.
#'  
#' @param dht.results a list of objects of class dht
#' @param species.name a character vectors detailing the 
#'   species codes 
#' @param clusters boolean whether observations are clusters of individuals
#' @return a list of results with an element for each species
#' @note Internal function not intended to be called by user.
#' @author Laura Marshall
#'
format.dht.results <- function(dht.results, species.name, clusters){
# format.dht.results function to fromat the dht results
#
# Arguments:
#   dht.results      - a list of objects of class dht
#   species.name     - a character vector containing the species names
#
# Value:
#   a list of proprated results with an element for each species
#
# Function Calls: none
#
  #get all identified 
  identified.codes <- species.name
  
  #create results object to store results in
  results <- list()
  #add all data for identified codes to results
  for(id in seq(along = identified.codes)){
    results[[identified.codes[id]]]$individual$summary <- dht.results[[identified.codes[id]]]$individuals$summary[,c("Region", "Area", "CoveredArea", "Effort", "n", "ER")]
    results[[identified.codes[id]]]$individual$N <- cbind(dht.results[[identified.codes[id]]]$individuals$N[,c("Label","Estimate","df")], PercentUnidentified = rep(0,nrow(dht.results[[identified.codes[id]]]$individuals$N)))
    if(clusters){
      results[[identified.codes[id]]]$clusters$summary <- dht.results[[identified.codes[id]]]$clusters$summary[,c("Region", "Area", "CoveredArea", "Effort", "n", "k", "ER")]
      results[[identified.codes[id]]]$clusters$N <- cbind(dht.results[[identified.codes[id]]]$clusters$N[,c("Label","Estimate","df")], PercentUnidentified = rep(0,nrow(dht.results[[identified.codes[id]]]$clusters$N)))
      results[[identified.codes[id]]]$Expected.S <- dht.results[[identified.codes[id]]]$Expected.S[,c("Region", "Expected.S")]
    }
  }
  #add in new cluster size col same as existing cluster size col
  if(clusters){
    for(r in seq(along = results)){
      results[[identified.codes[r]]]$Expected.S$new.Expected.S <- results[[identified.codes[r]]]$Expected.S$Expected.S 
    }   
  }
  #return results        
  return(results)
}



