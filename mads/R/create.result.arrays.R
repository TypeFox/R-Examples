#' Creates a list of arrays for storing the dht results
#'
#' Creates a list of arrays. These are used to store the summary, abundance
#' and density outputs of the \code{dht} routine called from \code{mrds}.
#'  
#' @param species.name a list of all the species in the analysis
#' @param species.code.definitions a list with an element for each 
#'   unidentified code which contains a vector of corresponding identified 
#'   species codes or NULL if not required 
#' @param region.table dataframe of region records - Region.Label and Area
#' @param clusters boolean, TRUE if observations are of cluster, FALSE if
#'   observations are of individuals.
#' @param n the number of bootstrap iterations to be completed.
#' @return list of arrays
#' @note Internal function not intended to be called by user.
#' @author Laura Marshall
#' @keywords utility
#'
create.result.arrays <- function(species.name, species.code.definitions, region.table, clusters, n){

  identified.species <- NULL
  for(sp in seq(along = species.name)){
    if(length(species.code.definitions[[species.name[sp]]]) == 1){
      identified.species <- c(identified.species, species.name[sp]) 
    }
  }
  no.id.species   <- length(identified.species)                                       
  strata.name  <- as.character(region.table$Region.Label)
  no.strata    <- length(strata.name)
  if(no.strata == 1){
    strata.name <- NULL
    no.strata <- 0
  }
  #create arrays to record bootstrap results
  individual.summary <- array(dim=c(no.strata+1, 5, n, no.id.species), dimnames = list(c(strata.name, "Total"), c("Area", "CoveredArea", "Effort", "n", "ER"), 1:n, identified.species))                                    
  individual.N       <- array(dim=c(no.strata+1, 3, n, no.id.species), dimnames = list(c(strata.name, "Total"), c("Estimate", "df", "PercentUnidentified"), 1:n, identified.species))        
  if(!clusters){
    # NO CLUSTERS - store these arrays in a list
    bootstrap.results <- list(individual.summary = individual.summary, individual.N = individual.N) 
  }else{ 
    # CLUSTERS
    clusters.summary   <- array(dim=c(no.strata+1, 6, n, no.id.species), dimnames = list(c(strata.name, "Total"), c("Area", "CoveredArea", "Effort", "n", "k", "ER"), 1:n, identified.species))
    clusters.N         <- array(dim=c(no.strata+1, 3, n, no.id.species), dimnames = list(c(strata.name, "Total"), c("Estimate", "df", "PercentUnidentified"), 1:n, identified.species))
    Expected.S         <- array(dim=c(no.strata+1, 2, n, no.id.species), dimnames = list(c(strata.name, "Total"), c("Expected.S", "new.Expected.S"), 1:n, identified.species))
    # store these arrays in a list  
    bootstrap.results <- list(individual.summary = individual.summary, individual.N = individual.N, clusters.summary = clusters.summary, clusters.N = clusters.N, Expected.S = Expected.S)
  } 
  return(bootstrap.results)                               
}