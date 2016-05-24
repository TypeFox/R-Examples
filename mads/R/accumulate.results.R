#' Enters the prorated results into the bootstrap.results array 
#'  
#' @param n index of the current bootstrap iteration
#' @param bootstrap.results list of 4-dimensional arrays containing the 
#'   bootstrap results
#' @param formatted.results list of data objects similar to the dht class
#' @param clusters boolean are the observations clusters of individuals
#'   bootstrap results 
#' @return list of 4-dimensional arrays containing the updated 
#' @author Laura Marshall
#'
accumulate.results <- function(n, bootstrap.results, formatted.results, clusters){
# 
# execute.multi.analysis  - function for dealing with model uncertainty, covariate uncertainty and unidentified species in Distance Sampling
#
# Arguments:
#   n                  - iteration number
#   bootstrap.results  - list of 4-dimensional arrays
#   prorated.results   - list of data objects similar to the dht class 
#
# Value:
#   returns the updated list of 4-dimensional arrays
#
# Function Calls: none
#                                      
  species.name <- names(formatted.results)
  #fill in results for each species
  for(sp in seq(along=species.name)){
    #Fill in individual results
    bootstrap.results$individual.summary[,,n,species.name[sp]] <- as.matrix(formatted.results[[species.name[sp]]]$individual$summary[,2:6])
    bootstrap.results$individual.N[,,n,species.name[sp]] <- as.matrix(formatted.results[[species.name[sp]]]$individual$N[,2:4])
  }
  if(clusters){
    for(sp in seq(along=species.name)){
      #Fill in individual results
      bootstrap.results$clusters.summary[,,n,species.name[sp]] <- as.matrix(formatted.results[[species.name[sp]]]$clusters$summary[,2:7])
      bootstrap.results$clusters.N[,,n,species.name[sp]] <- as.matrix(formatted.results[[species.name[sp]]]$clusters$N[,2:4])
      #this can be changed back after dht is fixed
      #bootstrap.results$Expected.S[,,n,species.name[sp]] <- as.matrix(formatted.results[[species.name[sp]]]$Expected.S[,2:3])
      bootstrap.results$Expected.S[,,n,species.name[sp]] <- as.matrix(formatted.results[[species.name[sp]]]$Expected.S[1:(dim(as.matrix(formatted.results[[species.name[sp]]]$clusters$N[,2:4]))[1]),2:3])
    }
  }
  return(bootstrap.results)
}




