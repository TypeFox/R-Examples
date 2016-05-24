#
# create.population.R
#
# Copyright (c) 2010-2013 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified May, 2013
# first written Mar, 2011
# Contains: create.population
#

#  create.population
#
# DESCRIPTION:
#  Creating an object of class population using data supplied by user
# PARAMETERS:
#   - offspring$phenotypes - matrix containing offspring phenotype data (have to be supported, if not - function
#   - quits with error
#   - founders - matrix containing founders phenotype data (optional)
#   - offspring$genotypes - matrix containing offspring genotype data (optional)
#   - maps$genetic - matrix containing genetic map (optional)
#   - maps$physical - matrix containing physical map (optional)
# OUTPUT:
#  An object of class population
#

create.population <- function(offspringPhenotypes, founders, foundersGroups, offspringGenotypes, mapsGenetic, mapsPhysical,
  populationType=c("riself", "f2", "bc", "risib"), noWarn=FALSE, verbose=FALSE, debugMode=0){
  if(verbose && debugMode==1) cat("create.population starting.\n")
  s <- proc.time()
  population <- NULL
  populationType <- match.arg(populationType)
  
  if(missing(offspringPhenotypes)) stop("No offspring phenotype data provided!")
  population <- add.to.populationSub.internal(population, populationType, offspringPhenotypes, "offspring$phenotypes", verbose, debugMode)
  
  if(missing(founders)){
    population <- simulateParentalPhenotypes(population, population$offspring$phenotypes, populationType)
  }else{
    n.childrenNotInParental <- sum(!(rownames(founders)%in%rownames(population$offspring$phenotypes)))
    if(n.childrenNotInParental == nrow(founders)) stop("No match between the row names in the founders and offspring.\n")
    if(n.childrenNotInParental != 0) warning(n.childrenNotInParental,"markers from founders file are not present in offspring data and will be removed.\n")
    founders <- founders[which((rownames(founders) %in% rownames(population$offspring$phenotypes))),]
    population <- add.to.populationSub.internal(population, populationType, founders, "founders", verbose, debugMode)
  }
  
  if(missing(foundersGroups)) stop("No information about founders groups provided!\n")
  if(length(foundersGroups) !=ncol(population$founders$phenotypes)) stop("foundersGroup parameter should have length equal to number of columns in founders phenotype data") 
  population <- add.to.populationSub.internal(population, populationType, foundersGroups, "founders$groups", verbose, debugMode)

  if(missing(offspringGenotypes)){
    if(verbose && !(noWarn))cat("No offspring genotypic data provided. You can supply it later using add.to.population.\n")
  }else{
    population <- add.to.populationSub.internal(population, populationType, offspringGenotypes, "offspring$genotypes", verbose, debugMode)
  }
  if(missing(mapsGenetic)){
    if(verbose && !(noWarn))cat("No genotic map provided. You can supply it later using add.to.population.\n")
  }else{
    population <- add.to.populationSub.internal(population, populationType, mapsGenetic, "maps$genetic", verbose, debugMode)
  }
  if(missing(mapsPhysical)){
    if(verbose && !(noWarn))cat("No physical map provided.  You can supply it later using add.to.population.\n")
  }else{
    population <- add.to.populationSub.internal(population, populationType, mapsPhysical, "maps$physical", verbose, debugMode)
  }
  
  if(is.null(population)) stop("No data provided")
  class(population) <- c("population", populationType)
  check.population(population)
  if(verbose){
    e <- proc.time()
    cat("\ncreate.population finished")
    if(debugMode==2) cat(" in:",(e-s)[3],"seconds.\n")
  }
  invisible(population)
}