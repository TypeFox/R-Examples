##############################################################################
#                                                                            #
#                        GENETIC ALGORITHMS in R                             #
#                                                                            #      
#       This functions are the genectic algorithms used by the EFS           #
#             algorithms that are available in the package.                  #
#                                                                            #
##############################################################################

#
# 
# Genetic algorithm for SDIGA
#
#

.gaSDIGA <- function(type = c("binary", "real-valued", "permutation"), 
               fitness, ...,
               min, max, nBits,
               population ,
               selection,
               crossover, 
               mutation,
               popSize = 50, 
               pcrossover = 0.8, 
               pmutation = 0.1, 
               elitism = base::max(1, round(popSize*0.05)), 
               maxiter = 100,
               run = maxiter,
               maxfitness = Inf,
               names = NULL,
               suggestions = NULL, 
               keepBest = FALSE,
               parallel = FALSE,
               monitor = NULL,
               DNFRules = FALSE,
               seed = NULL) 
{
  
  call <- match.call()
  
  type <- match.arg(type)

  
  if(missing(fitness))
  { stop("A fitness function must be provided") }
  if(!is.function(fitness)) 
  { stop("A fitness function must be provided") }
  if(popSize < 10) 
  { warning("The population size is less than 10.") }
  if(maxiter < 1) 
  { stop("The maximum number of iterations must be at least 1.") }
  if(elitism > popSize) 
  { stop("The elitism cannot be larger that population size.") }
  if(pcrossover < 0 | pcrossover > 1)
  { stop("Probability of crossover must be between 0 and 1.") }
  if(is.numeric(pmutation))
  { if(pmutation < 0 | pmutation > 1)
  { stop("If numeric probability of mutation must be between 0 and 1.") }
  else if(!is.function(population))
  { stop("pmutation must be a numeric value in (0,1) or a function.") }
  }
  if(missing(min) & missing(max) & missing(nBits))
  { stop("A min and max range of values (for 'real-valued' or 'permutation' GA) or nBits (for 'binary' GA) must be provided!") }
  
  switch(type, 
         "binary"      = { nBits <- as.vector(nBits)[1]
                           #min <- max <- NA
                           nvars <- nBits 
         },
         "real-valued" = { min <- as.vector(min)
                           max <- as.vector(max)
                           nBits <- NA
                           if(length(min) != length(max))
                           { stop("min and max must be vector of the same length!") }
                           nvars <- length(max) 
         },
         "permutation" = { min <- as.vector(min)[1]
                           max <- as.vector(max)[1]
                           nBits <- NA
                           nvars <- length(seq(min,max)) 
         }
  )
  
 
  

  Fitness <- rep(NA, popSize + 2)
  
  
  
  object <- new(".ga", 
                call = call, 
                type = type,
                min = min, 
                max = max, 
                nBits = nBits, 
                names = if(is.null(names)) character() else names,
                popSize = popSize,
                iter = 0, 
                run = 1, 
                maxiter = maxiter,
                suggestions = matrix(),
                population = matrix(), 
                elitism = elitism, 
                pcrossover = pcrossover, 
                pmutation = if(is.numeric(pmutation)) pmutation else NA,
                fitness = Fitness, 
                summary = matrix(),
                bestSol = list()
                )
  n_evals <- 0
  if(!is.null(seed)) set.seed(seed)
  
  if(!DNFRules)
    nGenes <- nvars * popSize
  else
    nGenes <- (length(max) - 1) * popSize
  
  numMutaciones <- ceiling(pmutation * nGenes)
  
  object@DNFRules <- DNFRules
  object@maxValuesRule <- max 
  object@popNew <- matrix(as.double(NA), nrow = popSize, ncol = nvars)
  
  # generate beginning population
 
  Pop <- matrix(as.double(NA), nrow = popSize + 2, ncol = nvars)
  
 
  
   if( ! DNFRules) 
    Pop[1:popSize,] <- population(object)[1:popSize,]
    else 
    Pop[1:popSize,] <- population(object)[1:popSize,]
  
  object@population <- Pop
  
  # start iterations
  for(iter in seq_len(maxiter))
  {
    # evalute fitness function (when needed) 
   for(i in seq_len(popSize + 2))
      if(is.na(Fitness[i]))
      { Fitness[i] <- fitness(Pop[i,], ...) 
        n_evals <- n_evals + 1} 
   



# update object
object@iter <- iter
object@population <- Pop
object@fitness <- Fitness



ord <- order(Fitness, decreasing = TRUE)
PopSorted <- Pop[ord,,drop=FALSE]
FitnessSorted <- Fitness[ord]

#Keep the population sorted by fitness
object@population <- PopSorted
object@fitness <- FitnessSorted


# check stopping criteria

if(n_evals >= run) break  
if(max(Fitness, na.rm = TRUE) >= maxfitness) break
if(object@iter == maxiter) break  


# selection
sel <- selection(object)
PopSorted <- sel$population
FitnessSorted <- sel$fitness

  
object@population <- PopSorted
object@fitness <- FitnessSorted


# crossover Only cross the 2 best individuals

parents <- c(1,2)
Crossover <- crossover(object, parents) # Only the best individuals are crossed
PopSorted[popSize + parents,] <- Crossover$children
FitnessSorted[popSize + parents] <- Crossover$fitness

# mutation (only .mutate popLength * probMut chromosomes)
pm <- if(is.function(pmutation)) pmutation(object) else pmutation
if(DNFRules) nvars <- length(max) - 1
if(is.function(mutation) & pm > 0)
{ 
  genes <- sample(x = seq_len(nGenes), size = numMutaciones, replace = TRUE)
  cromosomas <- floor(genes / nvars) + 1
  vars <- (cromosomas %% nvars) + 1
  if(!DNFRules)
    Mutation <- matrix(data = NA, nrow = numMutaciones, ncol = nvars)
  else
    Mutation <- matrix(data = NA, nrow = numMutaciones, ncol = nBits)
  FitnessSorted[cromosomas] <- NA
  for(i in seq_len(length(vars))) 
  {     
    Mutation[i,] <- mutation(object, cromosomas[i], vars[i])
  }
  PopSorted[cromosomas,] <- Mutation 
  
  object@population <- PopSorted
  object@fitness <- FitnessSorted
}

Pop <- PopSorted
Fitness <- FitnessSorted

}


# get solution(s)
object@fitnessValue <- max(object@fitness, na.rm = TRUE)

# return an object of class 'ga'
return(object)
}





















#
#
# Genetic algorithm for MESDIF
#
#

.gaMESDIF <- function(type = c("binary", "real-valued", "permutation"), 
               fitness, ...,
               min, max, nBits,
               population ,
               selection ,
               crossover , 
               mutation ,
               popSize = 50, 
               pcrossover = 0.8, 
               pmutation = 0.1, 
               elitism = base::max(1, round(popSize*0.05)), 
               maxiter = 100,
               run = maxiter,
               maxfitness = Inf,
               names = NULL,
               suggestions = NULL, 
               keepBest = FALSE,
               parallel = FALSE,
               monitor = NULL,
               DNFRules = FALSE,
               seed = NULL) 
{
  
  call <- match.call()
  
  type <- match.arg(type)

  
  if(missing(fitness))
  { stop("A fitness function must be provided") }
  if(!is.function(fitness)) 
  { stop("A fitness function must be provided") }
  if(popSize < 10) 
  { warning("The population size is less than 10.") }
  if(maxiter < 1) 
  { stop("The maximum number of iterations must be at least 1.") }
  if(elitism > popSize) 
  { stop("The elitism cannot be larger that population size.") }
  if(pcrossover < 0 | pcrossover > 1)
  { stop("Probability of crossover must be between 0 and 1.") }
  if(is.numeric(pmutation))
  { if(pmutation < 0 | pmutation > 1)
  { stop("If numeric probability of mutation must be between 0 and 1.") }
  else if(!is.function(population))
  { stop("pmutation must be a numeric value in (0,1) or a function.") }
  }
  if(missing(min) & missing(max) & missing(nBits))
  { stop("A min and max range of values (for 'real-valued' or 'permutation' GA) or nBits (for 'binary' GA) must be provided!") }
  
  switch(type, 
         "binary"      = { nBits <- as.vector(nBits)[1]
                           #min <- max <- NA
                           nvars <- nBits 
         },
         "real-valued" = { min <- as.vector(min)
                           max <- as.vector(max)
                           nBits <- NA
                           if(length(min) != length(max))
                           { stop("min and max must be vector of the same length!") }
                           nvars <- length(max) 
         },
         "permutation" = { min <- as.vector(min)[1]
                           max <- as.vector(max)[1]
                           nBits <- NA
                           nvars <- length(seq(min,max)) 
         }
  )
  
  
 
  #Define Fitness as a matrix, in MESDIF, each value of objectives is evaluated individually
  Fitness <- matrix(NA, nrow = popSize + elitism, ncol = 4)
  #This counts the number of indivuals that are dominated by this one
  Dominados <- numeric(popSize + elitism)
  
  
  object <- new(".ga", 
                call = call, 
                type = type,
                min = min, 
                max = max, 
                nBits = nBits, 
                names = if(is.null(names)) character() else names,
                popSize = popSize,
                iter = 0, 
                run = 1, 
                maxiter = maxiter,
                suggestions = matrix(),
                population = matrix(), 
                elitism = elitism, 
                pcrossover = pcrossover, 
                pmutation = if(is.numeric(pmutation)) pmutation else NA,
                fitness = Fitness, 
                summary = matrix(),
                bestSol = list())
  
  #This GA runs until a number of evaluations is reached, not iterations !
  n_evals <- 0
  if(!is.null(seed)) set.seed(seed)
  
  #Compute the number of genes, the mutation probability is applied over the gene.
  if(!DNFRules)
    nGenes <- nvars * popSize
  else
    nGenes <- (length(max) - 1) * popSize
  
  numMutaciones <- floor(pmutation * nGenes)
  
  object@DNFRules <- DNFRules
  object@maxValuesRule <- max 
  object@popNew <- matrix(as.double(NA), nrow = popSize, ncol = nvars)
  
  # generate beginning population and elite population
  Pop <- matrix(as.double(NA), nrow = popSize, ncol = nvars)
  elitePop <- matrix(as.double(NA), nrow = elitism, ncol = nvars)
  
  #Sets the next gene to mute
  Mu_next <- ceiling(log(runif(1)) / log(1 - pmutation))

  
  


  # fill the rest with a random population
   if( ! DNFRules) 
    Pop[1:popSize,] <- population(object, 0.25, round(nvars*0.25))[1:popSize,]
    else 
      Pop[1:popSize,] <- population(object, 0.25, round((length(max) - 1)*0.25))[1:popSize,]

  object@population <- Pop
  
  NonDominated <- logical(popSize + elitism) #Indicate wheter an individual is non-dominated
  WhoDominateMe <- vector(mode = "list", length = popSize + elitism) #Indicate the inviduals which domain this one
 # AdaptationValue <- numeric(popSize + elitism) 
  
  dots <- list(...) #Catch dots arguments
  nObjs <- length( which(!is.na(dots[[9]])) ) - 1
  
  volumenEsfera <- .volSphere(nObjs)
  
  
  nvariables <- nvars
  # start iterations
  for(iter in seq_len(maxiter))
  {
    #Create the union of this two populations and reinitializate all values
    NonDominated[] <- FALSE
    Dominados[] <- 0
    
    UnionPop <- matrix(NA, nrow = popSize + elitism, ncol = nvariables)
   
    UnionPop[seq_len(NROW(Pop)), ] <- Pop
    UnionPop[(NROW(Pop) + 1):NROW(UnionPop),] <- elitePop
    
    #normalize DNF RULES
#     if(DNFRules){
#       UnionPop <- matrix(unlist(apply(X = UnionPop, MARGIN = 1, FUN = .normalizeDNFRule, max)), ncol = nBits, byrow = TRUE)
#     }
    
    #Remove duplicated individuals in the UnionPop
    UnionPop <- na.exclude(UnionPop)
    duplicados <- which(! duplicated(UnionPop))
    UnionPop <- UnionPop[duplicados, , drop = F]
    Fitness <- Fitness[duplicados, , drop = F]


    # evalute fitness function (when needed) 

    for(i in seq_len( NROW(UnionPop) ))
      if(all(is.na(Fitness[i,])))
      { Fitness[i,] <- fitness(UnionPop[i,], ...) 
        n_evals <- n_evals + 1} 
    
    

#Compute dominated and non-dominated rules and initial adaptation Value
#Numero de individuos a los que domina cada regla
f <- na.exclude(Fitness)[, seq_len(nObjs), drop = F]
n_Ind <- NROW(f)
for(i in seq_len(n_Ind)){
  nd <- apply(X = f, MARGIN = 1, FUN = function(x, regla){ all(regla <= x) & any(regla < x)}, f[i,])
  Dominados[i] <- sum( apply(X = f, MARGIN = 1, FUN = function(x, regla){ all(regla >= x) & any(regla > x)}, f[i,]) )
  NonDominated[i] <- all( ! nd )
  WhoDominateMe[[i]] <- which(nd)
}

#Initial Adaptation Value
AdaptationValue <- vapply(X = WhoDominateMe[seq_len(n_Ind)], FUN = function(x, dominados) sum(dominados[x]) ,  FUN.VALUE = 2, Dominados)
kthDistances <- numeric(length(AdaptationValue))
#Distance measurement (Distancia entre valores de la regla o entre valores de fitnes? Uso valores DE FITNESS)
distancia <- as.matrix( dist(x = f, method = "euclidean", upper = TRUE, diag = TRUE))^2

#Order distance (The first colum is the distance with respect himself !! )
# And calculate the final adaptation value for each rule
kTH <- floor( sqrt(n_Ind - 1) )
for(i in seq_len(NROW(distancia))){
  distancia[i,] <- distancia[i, order(distancia[i, ]), drop = F]
  
  #Gets k-th closest neighbor , in this case, k = sqrt(popSize - 1)
  if(distancia[i, kTH] == 0){
    #If k-th closest is 0, get the next closest greater than 0
    aux <- distancia[i, (kTH + 1):ncol(distancia)] > 0
    closest <- which(aux, useNames = F)
    if(length(closest) > 0)
      dist <- distancia[i, closest[1] + kTH]
    else 
      #Exception: All individuals are equal
      dist <- 1
  } else {
    dist <- distancia[i, kTH]
  }
  
  kthDistances[i] <-  1 / dist^nObjs * kTH / n_Ind / volumenEsfera
}

#Normalize kthDistances
kthDistances <- kthDistances / sum(kthDistances)
AdaptationValue <- AdaptationValue + kthDistances
  
  
  
#Fill elite population
#Compute the number of non-dominated individuals
cantidadNoDominados <- sum(NonDominated)

if( cantidadNoDominados <= elitism){
  
  #Adition operator (Incluye los elitism mejores valores de adaptacion  a la poblacion elite)
  eliteIndividuals <- order(AdaptationValue)[seq_len(elitism)]
  elitePop <- UnionPop[eliteIndividuals, , drop = F]
 
} else {
  #Truncation operator
  lista <- .truncOperator(NonDominatedPop = UnionPop[which(NonDominated), , drop = F], elitePopSize = elitism, FitnessND = Fitness[which(NonDominated), , drop =  F])
  elitePop <- lista[[1]]
  eliteIndividuals <- which(NonDominated)[lista[[2]]]
  
  } 

# update object
object@iter <- iter
object@population <- Pop
object@fitness <- Fitness



# check stopping criteria

if(n_evals >= run) break  
if(max(Fitness, na.rm = TRUE) >= maxfitness) break
if(object@iter == maxiter) break  


# selection by Binary Tournament (Adaptation Values is the value for the "fitness")
#Copy the selected populatin into an intermediary population
sel <- selection(elitePop, popSize, nvariables, AdaptationValue[eliteIndividuals], Fitness[eliteIndividuals, , drop = F])
interPop <- sel$population
AdaptationValue <- sel$fitness
Fitness <- sel$obj

object@population <- interPop


# crossover performed by a double-point crossover
      if(pcrossover > 0)
        { 
          nmating <- round( (popSize/2) * pcrossover )
          
          #Create the population where we allocate the descendents from crossover and mutation
          descPop <- matrix(NA, ncol = nvariables, nrow = popSize*2)
          
          #mating <- matrix(sample(seq_len(popSize), size = (2*nmating), replace = TRUE), ncol = 2, byrow = T)
          for(i in seq_len(nmating))
            { 
                parents <- sample(seq_len(popSize), size = 2, replace = TRUE)
                Crossover <- crossover(object, parents)
                descPop[c(2*i - 1, 2*i),] <- Crossover$children
              
                
            }             
          #object@population <- descPop
          #nGenes <- NROW(na.exclude(descPop)) * nvariables
        }

# mutation (only .mutate popLength * probMut chromosomes)
pm <- pmutation
if(DNFRules) nvars <- length(max) - 1
if(pm > 0)
{ 

suma <- nmating*2 + 1
while(Mu_next <= nGenes){
cromosoma <- ceiling( Mu_next  / nvars ) 
gen <- (Mu_next %% nvars) + 1

descPop[suma, ] <- mutation(object, cromosoma, gen)
suma <- suma + 1
#Calcuate next gene
Mu_next <- Mu_next + ceiling(log( runif(1) ) /  log(1 - pmutation))
}

Mu_next <- Mu_next - nGenes


#Replace the worst individuals in the population with the genereted in crossovers and mutations

orden <- order(AdaptationValue)
#orden <- .qsort(AdaptationValue, left = 1, right = length(AdaptationValue), index = seq_len(length(AdaptationValue)))
#orden <- orden$indices
Pop <- interPop#[orden, , drop = F]
orden <- c(orden, (popSize+1):NROW(Fitness))
#Fitness <- Fitness[orden,, drop = F] 


Pop[orden[popSize:(popSize - (suma - 2))], ] <- descPop[seq_len(suma - 1),]
Fitness[orden[popSize:(popSize - (suma - 2))], ] <- NA

  
  }
}

# Return Non-duplicated individuals in elite pop
  if(DNFRules){
    elitePop <- matrix(unlist(apply(X = elitePop, MARGIN = 1, FUN = .normalizeDNFRule, max)), ncol = nBits, byrow = TRUE)
    } 
    elitePop[which(!duplicated(elitePop)), , drop = F]
  
}
























#
#
# Genetic algorithm for NMEEF-SD
#
#

.gaNMEEF <- function(type = c("binary", "real-valued", "permutation"), 
                     fitness, ...,
                     min, max, nBits,
                     population ,
                     selection ,
                     crossover , 
                     mutation ,
                     popSize = 50, 
                     pcrossover = 0.8, 
                     pmutation = 0.1, 
                     elitism = base::max(1, round(popSize*0.05)), 
                     maxiter = 100,
                     run = maxiter,
                     maxfitness = Inf,
                     names = NULL,
                     suggestions = NULL, 
                     keepBest = FALSE,
                     parallel = FALSE,
                     monitor = NULL,
                     DNFRules = FALSE,
                     seed = NULL,
                     porcCob = 0.5,
                     StrictDominance = TRUE,
                     reInitPop = TRUE, 
                     minCnf = 0.6) 
{
  
  call <- match.call()
  
  type <- match.arg(type)

  
  if(missing(fitness))
  { stop("A fitness function must be provided") }
  if(!is.function(fitness)) 
  { stop("A fitness function must be provided") }
  if(popSize < 10) 
  { warning("The population size is less than 10.") }
  if(maxiter < 1) 
  { stop("The maximum number of iterations must be at least 1.") }
  if(elitism > popSize) 
  { stop("The elitism cannot be larger that population size.") }
  if(pcrossover < 0 | pcrossover > 1)
  { stop("Probability of crossover must be between 0 and 1.") }
  if(is.numeric(pmutation))
  { if(pmutation < 0 | pmutation > 1)
  { stop("If numeric probability of mutation must be between 0 and 1.") }
    else if(!is.function(population))
    { stop("pmutation must be a numeric value in (0,1) or a function.") }
  }
  if(missing(min) & missing(max) & missing(nBits))
  { stop("A min and max range of values (for 'real-valued' or 'permutation' GA) or nBits (for 'binary' GA) must be provided!") }
  
  switch(type, 
         "binary"      = { nBits <- as.vector(nBits)[1]
         #min <- max <- NA
         nvars <- nBits 
         },
         "real-valued" = { min <- as.vector(min)
         max <- as.vector(max)
         nBits <- NA
         if(length(min) != length(max))
         { stop("min and max must be vector of the same length!") }
         nvars <- length(max) 
         },
         "permutation" = { min <- as.vector(min)[1]
         max <- as.vector(max)[1]
         nBits <- NA
         nvars <- length(seq(min,max)) 
         }
  )
  
  
  
  #Define Fitness as a matrix, in NMEEF, each value of objectives is evaluated individually
  Fitness <- matrix(NA, nrow = popSize * 2, ncol = 4)

  
  
  object <- new(".ga", 
                call = call, 
                type = type,
                min = min, 
                max = max, 
                nBits = nBits, 
                names = if(is.null(names)) character() else names,
                popSize = popSize,
                iter = 0, 
                run = 1, 
                maxiter = maxiter,
                suggestions = matrix(),
                population = matrix(), 
                elitism = elitism, 
                pcrossover = pcrossover, 
                pmutation = if(is.numeric(pmutation)) pmutation else NA,
                fitness = Fitness, 
                summary = matrix(),
                bestSol = list())
  #This GA runs until a number of EVALUATIONS is reached, not iterations !
  n_evals <- 0
  if(!is.null(seed)) set.seed(seed)
  
  
  #Compute the number of genes, the mutation probability is applied over the gene.
  if(!DNFRules)
    nGenes <- nvars * popSize
  else
    nGenes <- (length(max) - 1) * popSize
  
  
  #Set the next gene to mute
  numMutaciones <- round(pmutation * nGenes)
  
  object@DNFRules <- DNFRules
  object@maxValuesRule <- max 
  object@popNew <- matrix(as.double(NA), nrow = popSize, ncol = nvars)
  
  # generate beginning population and Offspring pop
  Pop <- matrix(as.integer(NA), nrow = popSize, ncol = nvars)
  OffspringPop <- matrix(as.integer(NA), nrow = popSize, ncol = nvars)
  UnionPop <- matrix(as.integer(NA), nrow = popSize *2, ncol = nvars)
  
 
  if( ! DNFRules) 
    Pop[1:popSize,] <- population(object, 0.25, round(nvars*0.25))#[1:(popSize-ng),]
  else 
    Pop[1:popSize,] <- population(object, 0.25, round((length(max) - 1)*0.25))#[1:(popSize,]
  
  object@population <- Pop
  
  #This counts the number of indivuals that domain this one
  Dominados <- numeric(popSize *2)
  WhoIDomain <- vector(mode = "list", length = popSize *2) #Indicate the inviduals which are dominated by this one
  rank <- numeric(popSize * 2) #Inidicate the rank of the individual
  CrowdingDistance <- numeric(popSize * 2) #Indicate the crowding distance of the individual
  dots <- list(...) #Catch dots arguments
  nObjs <- sum(! is.na(dots[[9]])) - 1 # Number of objetives we are using. -1 because the last value indicate te use of dnf or can representation for fitness calc.
  frentes <- vector(mode = "list", length = popSize * 2)
  fitnessFrentes <- vector(mode = "list", length = popSize * 2)
  coveredByIndividualFrentes <- vector(mode = "list", length = popSize * 2)
  coveredByIndividual <- matrix(FALSE, ncol = popSize * 2, nrow = length( dots[[1]][["data"]] ))
  cubiertoActual <- cubiertoAnterior <- logical(length( dots[[1]][["data"]] ))
  nIterEvolve <- 0 #Evaluation where the population evolved the last time
  fivePercent <- floor(run * 0.05)
  nvariables <- nvars
  dataset <- matrix(unlist(dots[[1]][["data"]]), nrow = nvariables + 1)
  targetClass <- which(dots[[1]][["class_names"]] == dots[[3]]) - 1
  
  
 
  #Evaluation of pop
  for(i in seq_len( NROW(Pop) ))
    if(all(is.na(Fitness[i,]))){
      fit <- fitness(Pop[i,], ...) 
      Fitness[i,] <- fit[[1]]
      coveredByIndividual[,i] <- as.logical(fit[[2]])
      n_evals <- n_evals + 1
    } 
  
  
  # start iterations
  for(iter in seq_len(maxiter))
  {
    #reinitialize all values


    UnionPop[] <- NA
    OffspringPop[] <- NA
    
    
    
    
    #Check stopping criteria
    if(n_evals >= run) break  
    if(max(Fitness, na.rm = TRUE) >= maxfitness) break
    if(object@iter == maxiter) break  
    
    frentes <- vector(mode = "list", length = popSize * 2)
    
    # selection by Binary Tournament 
    # Copy the selected populatin into the offspring population
 
    sel <- selection(Pop, popSize, rank, CrowdingDistance, Fitness, coveredByIndividual) 
    OffspringPop <- sel[[1]]
    FitnessOffspring <- sel[[2]]
    coveredByIndividual[,(popSize + 1):(popSize*2)] <- sel[[3]]
   
    
    object@population <- OffspringPop
   
    # crossover performed by a double-point crossover on Offspring Pop
    if(pcrossover > 0)
    { 
      nmating <- round( (popSize/2) * pcrossover )
      
      #mating <- matrix(sample(seq_len(popSize), size = 2*popSize, replace = TRUE), ncol = 2)
      mating <- matrix(sample(seq_len(popSize), size = floor(popSize/2) * 2, replace = TRUE), ncol = 2)
      equals <- which(mating[,1] == mating[,2])
      
      while(length(equals) > 0){
        mating[equals,] <- matrix(sample(seq_len(popSize), size = 2*length(equals), replace = TRUE), ncol = 2)
        equals <- which(mating[,1] == mating[,2])
      }
      
      #throw popSize/2 random numbers
      dados <- runif(floor(popSize/2))
    
      
      #Check which pair of individuals cross
      mating <- mating[which(dados <= pcrossover), , drop = F]
    
      for(i in seq_len(NROW(mating)))
      { 
        parents <- mating[i,]
        Crossover <- crossover(object, parents)
        OffspringPop[parents,] <- Crossover$children
        FitnessOffspring[parents,] <- NA
        
      }             
    
      object@population <- OffspringPop
      nGenes <- NROW(na.exclude(OffspringPop)) * nvariables
      
    }
    
    # mutation 
    pm <- pmutation
    if(DNFRules) nvars <- length(max) - 1
    if(pm > 0)
    { 
      dados <- runif(nGenes)
   
      
      genes <- which(dados <= pmutation)
     
      cromosomas <- ceiling(genes / nvariables)
      vars <- (genes %% nvariables) + 1
      
      for(i in seq_len(length(vars))) 
      {     
        object@population[cromosomas[i],] <- mutation(object, cromosomas[i], vars[i])
      }
    
      OffspringPop <- object@population
      
      FitnessOffspring[cromosomas,] <- NA
    
      
    }
    
    #Copy FitnessOffspring and evaluate individuals crossed and mutated
    Fitness[(NROW(Pop) + 1):(popSize*2), ] <- FitnessOffspring
  
    
   #Generate the next population.
    
    #Combine Pop and  OffspringPop into UnionPop
    
    UnionPop[seq_len(NROW(Pop)), ] <- Pop
    UnionPop[(NROW(Pop) + 1):NROW(UnionPop), ] <- OffspringPop
    #Evaluation of pop
    for(i in (NROW(Pop) + 1):(popSize*2) )
      if(all(is.na(Fitness[i,]))){
        fit <- fitness(UnionPop[i,], ...) 
        Fitness[i,] <- fit[[1]]
        coveredByIndividual[,i] <- as.logical( fit[[2]] )
        n_evals <- n_evals + 1
      } 
   
  #Compute dominance values for performing fast sorting algorithm
    f <- na.exclude(Fitness)[,seq_len(nObjs),drop=F]
    n_Ind <- NROW(f)
    for(i in seq_len(n_Ind)){
      nd <- apply(X = f, MARGIN = 1, FUN = .calculateDominance, f[i,,drop=F], StrictDominance) #Sangria de tiempo
      Dominados[i] <- length(which(nd == 1L))
      WhoIDomain[[i]] <- which(nd <= 0L)
    }
    
    #Get the Pareto front
    frentes[[1]] <- which(Dominados == 0)
    
    #Order UnionPop by dominance fronts (fast sorting algorithm)
    p <- 1
    while(length(frentes[[p]]) != 0){
      p <- p+1
      
      for(i in frentes[[p - 1]]){
        
       if(length(WhoIDomain[[i]] > 0)){
        Dominados[ WhoIDomain[[i]] ] <- Dominados[ WhoIDomain[[i]] ] - 1
        alSaco <-  which(Dominados == 0)
        if(length(alSaco) > 0){
          frentes[[p]] <- c(frentes[[p]], alSaco)
          rank[alSaco] <- p-1
        }
       } 
      } 

    }

  
    #Check if non-dominated front covers new examples and evolve
    enElFrente <- frentes[[1]]
    cubiertoActual <- apply(X = coveredByIndividual[, enElFrente, drop = F], MARGIN = 1, FUN = any)
    
    #Check if there are new examples covered 
    evolve <- any(! cubiertoAnterior[which(cubiertoActual)])
    

    
    for(i in seq_len(p - 1)){
      fitnessFrentes[[i]] <- Fitness[frentes[[i]], ,drop = F]
      coveredByIndividualFrentes[[i]] <- coveredByIndividual[,frentes[[i]], drop = F]
      frentes[[i]] <- UnionPop[frentes[[i]], , drop = F]
      
    }

      if(evolve){
        nIterEvolve <- n_evals
      } 
    cubiertoAnterior <- cubiertoActual
    
      # fill the next population 
      aux <- .fillPopulation(frentes, p - 1, fitnessFrentes, coveredByIndividualFrentes, popSize, nObjs)
      Pop <- aux[[1]]
      CrowdingDistance[seq_len(popSize)] <- aux[[2]]
      Fitness[seq_len(popSize), ] <- aux[[3]]
      rank <- aux[[4]]
      coveredByIndividual <- aux[[5]]
      
      # Checks if we need to reinitialize the population
      if(! evolve & ! aux[[6]] & reInitPop){
      # Check reinit condicion
      if(n_evals - nIterEvolve >= fivePercent){
        
        #re-initialize population
    
        #sel <- .reInitPob(elitePop = frentes[[1]], fitnessElite = fitnessFrentes[[1]], coveredElite = coveredByIndividual[,enElFrente], calculateDominance = CrowdingDistance, pctVariables = 0.5, cubiertoActual = cubiertoActual, dataset = dataset, maxRegla = dots[[1]][["conjuntos"]], cate = dots[[14]], num = dots[[15]], crispSets = dots[[1]][["crispSets"]], targetClass = targetClass, popSize = popSize )
        sel <- .reInitPob(elitePop = Pop, fitnessElite = Fitness[seq_len(popSize), ], coveredElite = coveredByIndividual[, seq_len(popSize)], crowdingDistance = CrowdingDistance, pctVariables = porcCob, cubiertoActual = cubiertoActual, dataset = dataset, maxRegla = dots[[1]][["conjuntos"]], cate = dots[[14]], num = dots[[15]], crispSets = dots[[1]][["crispSets"]], targetClass = targetClass, popSize = popSize )
        Pop <- sel[[1]]
        Fitness <- sel[[2]]
        CrowdingDistance <- sel[[3]]
        coveredByIndividual[,seq_len(popSize) ] <- sel[[4]]
        
        #Evaluation of generated pop
        for(i in seq_len(popSize) )
          if(all(is.na(Fitness[i,]))){
            fit <- fitness(Pop[i,], ...) 
            Fitness[i,] <- fit[[1]]
            coveredByIndividual[,i] <- as.logical( fit[[2]] )
            newCovered <- which(! cubiertoAnterior[which(as.logical( fit[[2]] ))]) # Check if this rule covers new uncovered examples
            n_evals <- n_evals + 1
            if(length(newCovered) > 0){
              nIterEvolve <- n_evals
              cubiertoAnterior[newCovered] <- TRUE
            }
          }
        
        #Check if new population evolve
        cubiertoActual <- apply(X = coveredByIndividual[, enElFrente, drop = F], MARGIN = 1, FUN = any)
        evolve <- any(! cubiertoAnterior[which(cubiertoActual)])
        cubiertoAnterior <- cubiertoActual
        
        if(evolve)
          nIterEvolve <- n_evals
         
        rank[] <- 0
        
      }
      }
    
      object@population <- Pop
      object@fitness <- Fitness
    
  }
  
  
  
  
  #get the last ranking
  
  #Compute dominance values
  f <- na.exclude(Fitness)[1:popSize,seq_len(nObjs), drop = F]
  n_Ind <- NROW(f)
  for(i in seq_len(n_Ind)){
    nd <- apply(X = f, MARGIN = 1, FUN = .calculateDominance, f[i,], TRUE)
    Dominados[i] <- length(which(nd == 1L))
    WhoIDomain[[i]] <- which(nd <= 0L)
  }
  
  #Get the Pareto front
  frentes[[1]] <- which(Dominados == 0)
  frentes[[1]] <- Pop[frentes[[1]], , drop = F]

  
  #Return individuals of the Pareto that has more confidence than the minimum.
  unicos <- which(!duplicated(frentes[[1]]))
  frentes[[1]] <- frentes[[1]][unicos, , drop = F]
  
  #Evaluate indivuduals for getting fuzzy confidence
  dots[[9]] <- list(.confianzaDifusa, NA, NA, FALSE)
  
  for(i in seq_len(NROW( frentes[[1]])) ){
      fit <- fitness(frentes[[1]][i,], dots[[1]],dots[[2]],dots[[3]],dots[[4]],dots[[5]],dots[[6]],dots[[7]],dots[[8]],dots[[9]],dots[[10]],dots[[11]],dots[[12]],dots[[13]],dots[[14]],dots[[15]],dots[[16]]) 
      Fitness[i,4] <- fit[[1]][1]
  }
  
  frentes[[1]] <- frentes[[1]][which(Fitness[seq_len(NROW(frentes[[1]])),4] > minCnf), , drop = F]
  frentes[[1]] #Return
  
}







#
#
#  FuGePSD Genetic Algorithm
#
#

.gaFuGePSD <- function(type,           # Type of execution (1 for One vs All, != 1 for normal execution)
                       dataset,        # keel object asociated to this genetic Algorithm (training file)
                       selection ,     # Selection function !
                       crossover ,     # Crossover function !
                       mutation ,      # mutation function
                       popSize = 50,      #size of the population
                       pcrossover = 0.8,  #Crossover Probability
                       pmutation = 0.1,   #Mutation Probability
                       pinsertion = 0.05,    #Insertion Probability
                       pdropping = 0.05,     #Dropping Probability
                       selectionSize = 2, #Tournament selection size
                       AllClass = TRUE,   #ALL_CLASS Attribute
                       T_norm = 1,        # T-norm used
                       ruleWeight = 0 ,    # Rule Weighting method to use
                       frm = 0,           # Fuzzy Reasoning Method to use
                       maxiter = 100,    #Max generations to run this genetic Algorithm.
                       weightsGlobalFitness = c(0.25, 0.25, 0.25, 0.25), #Weights Used in population Global Evaluation
                       seed = .randInt(0, 20000000)
                      )
{
  #First of all, we must check types of all attributes
  if(class(dataset) != "keel")
    stop("'dataset' must be a keel dataset object.")
  if(! is.function(selection))
    stop("'selection' must be function.")
  if(! is.function(crossover))
    stop("'crossover' must be function.")
  if(! is.function(mutation))
    stop("'mutation' must be function.")  
  if(popSize <= 0)
    stop("'popSize' must be greater than zero.")
 
  if(selectionSize < 2)
    stop("'selectionSize' must be greater than 2.")
  if(! is.logical(AllClass))
    stop("'AllClass' must be a logical value.")
  if(maxiter < 1)
    stop("'maxiter' must be greater than zero")

  if(length(weightsGlobalFitness) != 4)
    stop("length of 'weightsGlobalFitness' must be 4")
  
  suma <- sum(pcrossover, pmutation, pinsertion, pdropping)
  if(suma != 1 ){
    pcrossover <- pcrossover / suma
    pmutation <- pmutation / suma
    pinsertion <- pinsertion / suma
    pdropping <- pdropping / suma
  }
  
  #Once checked, the evolutive process can start
  if(type == 0){
    #Execution One Vs all
    
  } else {
    #Normal execution and Return
    executionPSD( clas = NULL,
                  dataset,        
                  selection ,    
                  crossover ,     
                  mutation ,      
                  popSize,      
                  pcrossover, 
                  pmutation ,  
                  pinsertion,    
                  pdropping ,     
                  selectionSize, 
                  AllClass,   
                  T_norm,       
                  ruleWeight,    
                  frm,           
                  maxiter,   
                  weightsGlobalFitness,
                  seed)
  }
  
}



executionPSD <- function(clas = NULL,   # number of the class to generate rules.
                       dataset,        # keel object asociated to this genetic Algorithm (training file)
                       selection ,     # Selection function !
                       crossover ,     # Crossover function !
                       mutation ,      # mutation function
                       popSize = 50,      #size of the population
                       pcrossover = 0.8,  #Crossover Probability
                       pmutation = 0.1,   #Mutation Probability
                       pinsertion = 0.05,    #Insertion Probability
                       pdropping = 0.05,     #Dropping Probability
                       selectionSize = 2, #Tournament selection size
                       AllClass = TRUE,   #ALL_CLASS Attribute
                       T_norm = 1,        # T-norm used
                       ruleWeight = 0 ,    # Rule Weighting method to use
                       frm = 0,           # Fuzzy Reasoning Method to use
                       maxiter = 100,    #Max generations to run this genetic Algorithm.
                       weightsGlobalFitness = c(0.25, 0.25, 0.25, 0.25), #Weights Used in population Global Evaluation
                       seed = .randInt(0, 20000000)
){
  populationFitness <- bestPopulationFitness <- numeric(1)
  exampleClass <- unlist(.getClassAttributes(dataset$data))
  
  #get categorical and numerical variables
  categorical <- dataset$atributeTypes == "c"
  categorical <- categorical[-length(categorical)]
  numerical <- !categorical
  
  datasetNoClass <- matrix(unlist(.separar(dataset)), nrow = dataset$nVars, ncol = dataset$Ns)
  bestPop <- vector(mode = "list", length = popSize)
  
  #Init population
  pop <- lapply(seq_len(popSize), function(x, dataset, tnorm, tconorm, rule_weight, clase){
    createNewRule(dataset, tnorm, tconorm, rule_weight, clase)
  }, dataset, T_norm, T_norm, ruleWeight, clas)
  
  #evaluate initial population individuals (In parallel for Linux)
  if(length(pop) >= 20 & Sys.info()[1] == "Linux"){
    pop <- parallel::mclapply(pop, Rule.evaluate, dataset, datasetNoClass, categorical, numerical, T_norm, ruleWeight, mc.cores = parallel::detectCores())
  } else {
   pop <- lapply(pop, Rule.evaluate, dataset, datasetNoClass, categorical, numerical, T_norm, ruleWeight)
  }
  #evaluate the whole population
  populationFitness <- Pop.evaluate(pop, dataset, datasetNoClass, exampleClass, frm, categorical, numerical, T_norm, weightsGlobalFitness)
  
  #best population is now initial population.
  bestPop <- pop
  bestPopulationFitness <- populationFitness
  cat(paste("Global Fitness obtained in generation [0]:", bestPopulationFitness, "\n", sep = " "))
  
  
  
  #Init the evolutive process
  for(generation in seq_len(maxiter - 1)){
    #First, create a join population with twice length of population. Then add pop to joinPop
    joinPop <- vector(mode = "list", length = length(pop) * 2)
    joinPop[seq_len(length(pop))] <- pop
    
    #Now we need to generate an offspring population length equal to pop length.
    #This offspring population is generated via genetic operators.
    
    dados <- runif(length(pop))
    first_parents <- vapply(X = seq_len(length(pop)), 
                            FUN = function(x, pop, tam){
                                  tournamentSelection(pop, tam)}, numeric(1), pop, selectionSize)
    
    #Specify the genetic operator to apply according to their probability in 'dados'
    cruzan <- first_parents[which(dados < pcrossover)]
    mutan <- first_parents[which(pcrossover <= dados & dados < (pcrossover + pmutation))]
    insertan <- first_parents[which((pcrossover + pmutation) <= dados & dados < (pcrossover + pmutation + pinsertion))]
    dropean <- first_parents[which(pcrossover + pmutation + pinsertion <= dados)]
    
    posJoinPop <- length(pop) + 1
    #Make crossovers
    for(i in cruzan){
      second_parent <- .randIntExcluded(1, length(pop), i)
      joinPop[[posJoinPop]] <- FuGePSD_crossover(rule1 = pop[[i]], rule2 = pop[[second_parent]], nvars = dataset$nVars + 1)
      posJoinPop <- posJoinPop + 1  
    }
    
    #Make mutations 
    for(i in mutan){
      joinPop[[posJoinPop]] <- FuGePSD_Mutation(pop[[i]], dataset)
      posJoinPop <- posJoinPop + 1
    }
    
    #Make insertions
    for(i in insertan){
      if(length(pop[[i]][[1]]) == dataset[[6]]){ 
        #If we cannot add more variables, we introduce a rule with an empty antecedent.
        joinPop[[posJoinPop]] <- Rule.clearAntecedent(pop[[i]])
      } else {
        #Add a random variable
        joinPop[[posJoinPop]] <- Rule.addVariable(pop[[i]], dataset)
      }
      posJoinPop <- posJoinPop + 1
    }
    
    #Make droppings
    for(i in dropean){
      if(length(pop[[i]][[1]]) == 1){
        #We cannot delete more variables, return an empty rule
        joinPop[[posJoinPop]] <- Rule.clearAntecedent(pop[[i]])
      } else {
        joinPop[[posJoinPop]] <- Rule.deleteVariable(pop[[i]])
      }
      posJoinPop <- posJoinPop + 1
    }
    
    #Evaluate joinPop.
    joinPop <- lapply(joinPop, Rule.evaluate, dataset, datasetNoClass, categorical, numerical, T_norm, ruleWeight)
    
    #Apply Token Competition
    pop <- tokenCompetition(joinPop, dataset)
    
    #Evaluate Global Fitness
    populationFitness <- Pop.evaluate(pop, dataset, datasetNoClass, exampleClass, frm, categorical, numerical, T_norm, weightsGlobalFitness)
    
    #Substitute best population if actual if better.
    if(bestPopulationFitness < populationFitness){
      bestPopulationFitness <- populationFitness
      bestPop <- pop
      cat(paste("Global Fitness obtained in generation [", generation, "]: ", bestPopulationFitness, "\n", sep = ""))
    }
    
    #cat("\r", (generation / (maxiter-1)) * 100, "% Completed.", sep = "")
  }
  
  #Order bestPop by conf_f (desc. order)
  fuzzy_conf <- vapply(X = bestPop, FUN = function(x){x$qm_Cnf_f}, numeric(1))
  sens <- vapply(X = bestPop, FUN = function(x){x$qm_Sens}, numeric(1))
  orden <- order(fuzzy_conf, decreasing = TRUE)
  fuzzy_conf <- fuzzy_conf[orden]
  sens <- sens[orden]
  
  bestPop <- bestPop[orden]
  
  #Return 
  list(bestPop = bestPop, conf = fuzzy_conf, sensitivity = sens)

}


#-------------------------------------------------------------------------------
# ---  THis part is part of the definition of the "ga" class done in the GA Package
#--------------------------------------------------------------------------------
#
# We need to change this behaviour, we dont want to depend on this functions.
# Because we dont use it.



methods::setClassUnion(".numericOrNA", members = c("numeric", "logical", "matrix"))

methods::setClassUnion(".matrixOrList", members = c("matrix", "list"))

#Modification of the class 'ga' provided by the package "GA" created by Luca Scrucca.

methods::setClass(Class = ".ga", 
         representation(call = "language",
                        type = "character",
                        min = ".numericOrNA", 
                        max = ".numericOrNA", 
                        nBits = ".numericOrNA", 
                        names = "character",
                        popSize = "numeric",
                        iter = "numeric", 
                        run = "numeric", 
                        maxiter = "numeric",
                        suggestions = "matrix",
                        population = ".matrixOrList",
                        popNew = "matrix",
                        elitism = "numeric", 
                        pcrossover = "numeric", 
                        pmutation = ".numericOrNA",
                        fitness = ".numericOrNA",
                        summary = "matrix",
                        bestSol = "list",
                        fitnessValue = "numeric",
                        solution = "matrix",
                        maxValuesRule = ".numericOrNA",
                        DNFRules = "logical"
         ),
         package = "SDR" 
) 






#---------------------------------------------------------------------------------------------







## 
## Modification of Permutation ga operators.
##  This modification generate an integer random population in the range[min, max]
##

.generarPoblacion <- function(object, ...)
{
  # Generate a random permutation of size popSize in the range [min, max]  
  min <- object@min
  max <- object@max
  type <- object@type
  size <- object@popSize
  
  
  if(! object@DNFRules){ # real-valued indica que se usan reglas tipo CAN
    population <- matrix(as.double(NA), nrow = size, ncol = length( max ) )
    for(i in 1:size)
      for(j in 1:length(max))
        #Se genera la poblacion inicial, el valor de no participacion no cuenta ! 
        population[i,j] <- sample(0:(max[j]), size = 1, replace = TRUE)
    
  } else { # reglas DNF 
    v <- sample(x = c(0,1), size = max[length(max)] * size, replace = TRUE)
    population <- matrix(data = v, nrow = size, ncol = max[length(max)],byrow = TRUE)
  }
  
  return(population)
}






#  
# GENERA LA POBLACION INICIAL PARA MESDIF
# En el argumento ... deben de ir primero el porcentaje de poblacion que se genera completamente aleatorio
# y en segundo lugar el numero maximo de variables que participan en la regla.
.generarPoblacionMESDIF <- function(object, ...)
{
  # Generate a random permutation of size popSize in the range [min, max]  
  min <- object@min
  max <- object@max
  type <- object@type
  size <- object@popSize
  lista <- list(...)
  pctAleatorio <- lista[[1]]
  numVarMax <- lista[[2]]
  var_init <- logical(length(max))
  
  reglas <- ceiling(size * (1-pctAleatorio))
  if(! object@DNFRules){ # real-valued indica que se usan reglas tipo CAN
    population <- matrix(max, nrow = size, ncol = length( max ), byrow = TRUE )
    # Biased Init
    for(i in seq_len(reglas)){
      var_init[] <- F
      numVar <- sample(numVarMax, size = 1)
      for(j in seq_len(numVar)){
        var <- sample(length(max), size = 1)
        while(var_init[var]){  #Hay que incluir esto tambien a reglas DNF
          var <- sample(length(max), size = 1)
        }
        population[i, var] <- sample(0:(max[var]), size = 1, replace = TRUE) # No-Participate value is not into account
        var_init[var] <- T
        }
    }
    
    
    #Random Init
    for(i in (reglas + 1):size){
      
      for(j in seq_len(length(max)))
        
        population[i,j] <- sample(0:(max[j]), size = 1, replace = TRUE)
    }
    
    
  } else { # reglas DNF 
    
    # Random Init 
    v <- sample(x = c(0,1), size = max[length(max)] * size, replace = TRUE)
    population <- matrix(data = v, nrow = size, ncol = max[length(max)],byrow = TRUE)
    
    #Biased Init
    nReglasABorrar <- length(max) - 1 - numVarMax
    for(i in (reglas + 1):size){
      varia <- sample(length(max) - 1, size = nReglasABorrar, replace = FALSE)
      for(j in seq_len(nReglasABorrar)){
        population[i, ] <- .borrar_gen(regla = population[i,], variable = varia[j], max_valor_variables = max, DNF_Rules = TRUE)
      }
    }
  }
  return(population)
}









#
# Generates the initial population of NMEEF-SD
#

.generarPoblacionNMEEF <- function(object, ...)
{
  # Generate a random permutation of size popSize in the range [min, max]  
  min <- as.integer(object@min)
  max <- as.integer(object@max)
  type <- object@type
  size <- object@popSize
  lista <- list(...)
  pctAleatorio <- lista[[1]]
  numVarMax <- lista[[2]]
  var_init <- logical(length(max))
  
  reglas <- ceiling(size * (1-pctAleatorio))
  if(! object@DNFRules){ # real-valued indica que se usan reglas tipo CAN
    population <- matrix(max, nrow = size, ncol = length( max ), byrow = TRUE )
    # Biased Init
    for(i in seq_len(reglas)){
      var_init[] <- F
      numVar <- sample(numVarMax, size = 1)
      for(j in seq_len(numVar)){
        var <- sample(length(max), size = 1)
        while(var_init[var]){  #Hay que incluir esto tambien a reglas DNF
          var <- sample(length(max), size = 1)
        }
        population[i, var] <- sample(0:(max[var] - 1), size = 1, replace = TRUE) # No-Participate value is not into account
        var_init[var] <- T
      }
    }
      
      #Random Init
      for(i in (reglas + 1):size){
        
        for(j in seq_len(length(max)))
          
          population[i,j] <- sample(0:(max[j]), size = 1, replace = TRUE)
      }
      
  } else { # reglas DNF 
    
    # Random Init 
    v <- sample(x = c(0,1), size = max[length(max)] * size, replace = TRUE)
    population <- matrix(data = v, nrow = size, ncol = max[length(max)],byrow = TRUE)
    
    #Biased Init
    nReglasABorrar <- length(max) - 1 - numVarMax
    for(i in (reglas + 1):size){
      varia <- sample(length(max) - 1, size = nReglasABorrar, replace = FALSE)
      for(j in seq_len(nReglasABorrar)){
        population[i, ] <- .borrar_gen(regla = population[i,], variable = varia[j], max_valor_variables = max, DNF_Rules = TRUE)
      }
    }
  }

  return(population)
}









#
# Double-point crossover
#
.ga_dpCrossover <- function(object, parents, ...)
{


  if( ! object@DNFRules){ #REGLAS CAN
     
      parents <- object@population[parents,,drop = FALSE]
      n <- ncol(parents)
      children <- matrix(as.double(NA), nrow = 2, ncol = n)
      fitnessChildren <- rep(NA, 2)
      crossOverPoint1 <- sample(seq_len(n), size = 1, replace = TRUE)  
      if(crossOverPoint1 == (n) )
      { crossOverPoint2 <- n   } else {
      crossOverPoint2 <- sample((crossOverPoint1 + 1):n, size = 1, replace = TRUE)
      }
      
      children[1,] <- parents[1,]
      children[2,] <- parents[2,]
      
      
      
      children[1, crossOverPoint1:crossOverPoint2] <- parents[2, crossOverPoint1:crossOverPoint2]
      children[2, crossOverPoint1:crossOverPoint2] <- parents[1, crossOverPoint1:crossOverPoint2]
      
      out <- list(children = children, fitness = fitnessChildren)
      return(out)
  
      } else { # REGLAS DNF 
      
        
        parents <- object@population[parents,,drop = FALSE]
        max <- object@max
        n <- length(max)
        children <- matrix(as.double(NA), nrow = 2, ncol = max[length(max)])
        fitnessChildren <- rep(NA, 2)
        
        #Cambiar n - 1 por n, ya que es para comparar con java
        rangCrossover1 <- 2:(n) # Si la longitud de esto es 1, no se puede usar sample
        if(length(rangCrossover1) > 1){
          crossOverPoint1 <- sample(rangCrossover1, size = 1, replace = TRUE)  
        } else {
          crossOverPoint1 <- rangCrossover1
        }
        
        if(crossOverPoint1 == (n) )
        { crossOverPoint2 <- n   } else {
          rangCrossover2 <- (crossOverPoint1 + 1):n
          if(length(rangCrossover2) > 1){
            crossOverPoint2 <- sample(rangCrossover2, size = 1, replace = TRUE)
          } else {
            crossOverPoint2 <- rangCrossover2
          }
        }
      
          rango <- (max[crossOverPoint1 - 1] + 1):max[crossOverPoint2]
        
        
        children[1,] <- parents[1,]
        children[2,] <- parents[2,]
        
        children[1, rango] <- parents[2, rango]
        children[2, rango] <- parents[1, rango]
        
        out <- list(children = children, fitness = fitnessChildren)
        return(out)
        
    }
}









# SDIga Mutation operator

.gaCAN_Mutation <- function(object, parent, ...)
{
  
  mutar <- parent <- as.vector(object@population[parent,]) 
  mutar <-  .mutate(cromosoma = mutar, variable = ...[[1]], max_valor_variables = object@maxValuesRule, DNF_Rule = object@DNFRules )
  
  return(mutar)
}







# MESDIF Mutation operatos

..gaMESDIF_Mutation <- function(object, parent, ...)
{
  
  mutar <- parent <- as.vector(object@population[parent,]) 
  mutar <-  .mutateMESDIF(cromosoma = mutar, variable = ...[[1]], max_valor_variables = object@maxValuesRule, DNF_Rule = object@DNFRules )
  
  return(mutar)
}







#NMEEF-SD Mutation operator

..gaNMEEF_Mutation <- function(object, parent, ...)
{
  
  mutar <- parent <- as.vector(object@population[parent,]) 
  mutar <-  .mutateNMEEF(cromosoma = mutar, variable = ...[[1]], max_valor_variables = object@maxValuesRule, DNF_Rule = object@DNFRules )
  
  return(mutar)
}








#Selection Function of SDIGA

.ga_SDIgaSelection <- function(object){
  n <- object@popSize 
  object@population[(n + 1):nrow(object@population), ] <- NA
  
  
  object@fitness[(n + 1):length(object@fitness)] <- NA
  
  list(population = object@population, fitness = object@fitness)
  
}






# Binary Tournament selection operator for MESDIF

.ga_MESDIFBinTournamentSelection <- function(elitePop, sizePop, nvars, FitnessElite, ObjValues){
  newPop <- matrix(NA, nrow = sizePop, ncol = nvars)
 
  nas <- which(is.na(elitePop[,1,drop = F]))
  if(length(nas) > 0 )ObjValues <- ObjValues[- nas, , drop = F]
  elitePop <- na.exclude(elitePop)

  seleccion <- sample(NROW(elitePop), size = sizePop * 2, replace = TRUE)
  Fitness <- numeric(sizePop)
  Obj <- matrix(NA, nrow = sizePop + NROW(elitePop), ncol = 4)
  
  fit <- FitnessElite[seleccion] 
  sel <- matrix(seleccion, nrow = 2)
  fit <- matrix(fit, nrow = 2)
  
  vencedor <- fit[1,] <= fit[2,]

  num <- sum(vencedor)
 
  if(num > 0){
    b <- which(vencedor)
    a <- sel[1, b]
    newPop[b, ] <- elitePop[a, ]
    Fitness[b] <- FitnessElite[a]
    Obj[b, ] <- ObjValues[a, ]
  } 
  b <- which(!vencedor)
  a <- sel[2, b]
  newPop[b, ] <- elitePop[a, ]
  Fitness[b] <- FitnessElite[a]
  Obj[b, ] <- ObjValues[a, ]

  Obj[(sizePop + 1):NROW(Obj), ] <- ObjValues
  list(population = newPop, fitness = Fitness, obj = Obj)
  
  
}





#
#
# This function execute the corresponding genetic algorithm in function of the value of 'algorithm'
#
#

.ejecutarga <- function(algorithm, dataset, targetClass, n_vars, por_cubrir, nLabels, N_evals, tam_pob, p_cross = 0.5, p_mut, seed, Objetivos = c(.LocalSupport, .confianza, NULL, FALSE), Pesos = c(0.7,0.3,0), DNFRules = FALSE, cate, num, elitism = 5, porcCob = 0.5, strictDominance = TRUE, reInit = TRUE, minCnf = 0.6){
  
  ma <- dataset$conjuntos
  
if(DNFRules) {
  ma <- Reduce(f = '+', x = ma, accumulate = TRUE)
  ma <- c(0,ma)
}
  #Para reglas DNF, hay que utilizar como type el valor "binary", en vez de usar min y max, hay que usar el valor nBits, que indica la cantidad de genes que tiene cada cromosoma.
  #Tambien hay que utilizar el valor 'max' para saber cu?ntos genes pertenecen a cada variable.
  switch(algorithm, 
  "SDIGA" = { resultado <- .gaSDIGA(type = if(!DNFRules) "real-valued" else "binary", 
                  fitness = .fit13, dataset, matrix(unlist(.separar(dataset)), nrow = length(dataset[[2]]) - 1, ncol = length(dataset[[7]])), targetClass, por_cubrir, n_vars, nLabels, ma, FALSE, Objetivos, Pesos, DNFRules, Objetivos[[4]], FALSE,  cate, num,
                  min = dataset[[4]][-length(dataset[[4]])],
                  max = ma,
                  nBits = ma[length(ma)],
                  population = .generarPoblacion,
                  selection = .ga_SDIgaSelection,
                  crossover = .ga_dpCrossover,
                  mutation = .gaCAN_Mutation,  
                  popSize = tam_pob,
                  pcrossover = 1 / tam_pob, 
                  pmutation = p_mut, # / length(ma), #Mutation probability applied at the gene
                  elitism = 0,
                  maxiter = N_evals,#floor( (N_evals - tam_pob) / (2 + tam_pob  * p_mut)),
                  run = N_evals, # No queremos que se detenga la evaluacion.
                  #  maxfitness = 1, # Si encontramos un cromosoma que tiene valor maximo, detenemos la busqueda.
                  names = dataset[[2]][1:n_vars],
                  keepBest = FALSE,
                  parallel = FALSE,
                  monitor = NULL,
                  DNFRules = DNFRules,
                  seed = seed) 
  }, 
  "MESDIF" =  { resultado <- .gaMESDIF(type = if(!DNFRules) "real-valued" else "binary", 
                                #fitness = .fit12, dataset, .separar(dataset = dataset), targetClass, por_cubrir, n_vars, nLabels, ma, FALSE, Objetivos, Pesos, DNFRules, Objetivos[[4]], FALSE,  cate, num,# Parametros de .fit12
                                fitness = .fitnessMESDIF, dataset, matrix(unlist(.separar(dataset)), nrow = length(dataset[[2]]) - 1, ncol = length(dataset[[7]])), targetClass, por_cubrir, n_vars, nLabels, ma, FALSE, Objetivos, c(0.7, 0.3, 0), DNFRules, Objetivos[[4]], FALSE,  cate, num,
                                min = dataset[[4]][-length(dataset[[4]])],
                                max = ma,
                                nBits = ma[length(ma)],
                                population = .generarPoblacionMESDIF,
                                selection = .ga_MESDIFBinTournamentSelection,
                                crossover = .ga_dpCrossover,
                                mutation = ..gaMESDIF_Mutation,  
                                popSize = tam_pob,
                                pcrossover = p_cross, 
                                pmutation = p_mut / length(ma),
                                elitism = elitism,
                                maxiter = N_evals,
                                run = N_evals, # No queremos que se detenga la evaluacion.
                                names = dataset[[2]][1:n_vars],
                                keepBest = FALSE,
                                parallel = FALSE,
                                monitor = NULL,
                                DNFRules = DNFRules,
                                seed = seed) 
                return(resultado)
  }, 
  "NMEEFSD" =  { resultado <- .gaNMEEF(type = if(!DNFRules) "real-valued" else "binary", 
                                      #fitness = .fit12, dataset, .separar(dataset = dataset), targetClass, por_cubrir, n_vars, nLabels, ma, FALSE, Objetivos, Pesos, DNFRules, Objetivos[[4]], FALSE,  cate, num,# Parametros de .fit12
                                      fitness = .fitnessMESDIF, dataset, matrix(unlist(.separar(dataset)), nrow = length(dataset[[2]]) - 1, ncol = length(dataset[[7]])), targetClass, por_cubrir, n_vars, nLabels, ma, FALSE, Objetivos, c(0.7,0.3,0), DNFRules, Objetivos[[4]], FALSE,  cate, num, TRUE, 
                                      min = dataset[[4]][-length(dataset[[4]])],
                                      max = ma,
                                      nBits = ma[length(ma)],
                                      population = .generarPoblacionNMEEF,
                                      selection = .selectionNMEEF,
                                      crossover = .ga_dpCrossover,
                                      mutation = ..gaNMEEF_Mutation,  
                                      popSize = tam_pob,
                                      pcrossover = p_cross, 
                                      pmutation = p_mut,
                                      elitism = 0,
                                      maxiter = N_evals,
                                      run = N_evals, 
                                      names = dataset[[2]][1:n_vars],
                                      keepBest = FALSE,
                                      parallel = FALSE,
                                      monitor = NULL,
                                      DNFRules = DNFRules,
                                      seed = seed, 
                                      porcCob = porcCob,
                                      StrictDominance = strictDominance,
                                      reInitPop = reInit,
                                      minCnf = minCnf) 
  return(resultado) #Devolver las que superen minConf
  }
  )
  #Only for SDIGA
  .getBestRule(resultado)
}



#
# Mark examples of the dataset covered by the rule returned by genetic algorithm
# ONLY FOR SDIGA
#

.marcar_ejemplos <- function(regla, dataset, targetClass, nVars, maxRegla, por_cubrir, nLabels, Objetivos = c(.LocalSupport, .confianza, NA, FALSE), Pesos = c(0.7,0.3,0), DNFRules = FALSE, cate, num){
  
  
  #Devolver ejemplos nuevos cubiertos de la clase objetivo
  #Cambiar por fit13
  cover <- .fit13(regla = regla, dataset = dataset, noClass = matrix(unlist(.separar(dataset)), nrow = length(dataset[[2]]) - 1, ncol = length(dataset[[7]])), targetClass = targetClass, por_cubrir = por_cubrir, n_Vars = nVars,nLabels = nLabels, max_regla = maxRegla , marcar = TRUE, Objetivos = Objetivos, Pesos = Pesos, DNFRules = DNFRules, difuso = Objetivos[[4]], cate = cate, num = num)
  
  sumaNuevos <- sum(cover[[1]]) - sum(dataset[["covered"]])
  confi <- cover[[2]]
  por_cubrir <- por_cubrir - sumaNuevos
  
  return( list(cubreNuevos = sumaNuevos > 0, covered = cover , porCubrir = por_cubrir, confidence = confi) )
  
}












#
# Returns de best rule of the genetic algorithm. In case of draw, returns the rule with less atributes.
# 
# ONLY FOR SDIGA 
#
.getBestRule <- function(resultado){
  bestFitness <- resultado@fitnessValue
  
  empates <- which(resultado@fitness == bestFitness)
  if(length(empates) > 1){
    if(!resultado@DNFRules){
      lista <- apply(X = resultado@population[empates, ], MARGIN = 1, FUN = function(x, max) sum(x != max), resultado@max )
    } else{
      
      lista <- apply(X = resultado@population[empates, ], MARGIN = 1, FUN = function(x, max){
        particip <- .getParticipantes(regla = x, max_regla = max, DNFRules = TRUE)
        val <- numeric(0)
        for(i in 2:length(max)){
          if(particip[i-1])
            val <- c(val, (max[i-1]+1):max[i])
        }
        sum(x[val] != 0)
      } , resultado@max ) # ESTO HAY QUE CAMBIARLO !!
      
    }
    orden <- order(lista[which(lista > 0)])
    
    return(resultado@population[ orden[1] , ])
  } else {
    return(resultado@population[ empates , ])
  }
}















#Return the fitness value(s) of a rule

.fit13 <- function(regla, dataset, noClass, targetClass, por_cubrir, n_Vars, nLabels, max_regla, marcar = FALSE, Objetivos = c(.LocalSupport, .confianza, NULL, FALSE), Pesos = c(0.7,0.3,0), DNFRules = FALSE, difuso = FALSE ,test = FALSE, cate, num){

  
  if( ! any(is.na(regla))) { #Si la regla no tiene NA se puede evaluar
    
    regla <- as.integer(regla)
    participantes <- logical(length(max_regla))
    participantes <- .getParticipantes(regla = regla, max_regla = max_regla, DNFRules = DNFRules)
    
    
    #If it's not the empty rule
    if(any(participantes)){
      
      cat_particip <- which(cate & participantes)
      num_particip <- which(num & participantes)
      
      max_regla_cat <- max_regla[cat_particip]
      max_regla_num <- max_regla[num_particip]
      
      if(!DNFRules) { # CAN RULES
        
        #Split into numerical variables and categorical ones. (And participate in the rule)
        if(length(cat_particip) > 0){
          rule_cat <- regla[cat_particip]
        }
        
        if(length(num_particip) > 0){
          rule_num <- regla[num_particip]
          
          fuzzy_sets <- dataset[["fuzzySets"]][1:nLabels, 1:3, num_particip, drop = F]
          crispSets <- dataset[["crispSets"]][1:nLabels, 1:2, num_particip, drop = F]
          #  Get values for xmin, xmedio and xmax for fuzzy computation.   
          n_matrices <- dim(fuzzy_sets)[3]  
   
          xmin <- fuzzy_sets[cbind(rule_num + 1, 1, seq_len(n_matrices))]
          xmax <- fuzzy_sets[cbind(rule_num + 1, 3, seq_len(n_matrices))]
          xmedio <- fuzzy_sets[cbind(rule_num + 1, 2, seq_len(n_matrices))]
          
          #Get values for xmin and xmax for crisp computation
          n_matricesCrisp <- dim(crispSets)[3]  
          xminC <- crispSets[cbind(rule_num + 1, 1, seq_len(n_matricesCrisp))]
          xmaxC <- crispSets[cbind(rule_num + 1, 2, seq_len(n_matricesCrisp))]
        }
        
        gr_perts <- .compara_CAN9(ejemplo = noClass, rule_cat = rule_cat, rule_num = rule_num, catParticip = cat_particip, numParticip = num_particip, xmin = xmin, xmedio = xmedio, xmax = xmax, n_matrices = n_matrices, xminCrisp = xminC, xmaxCrisp = xmaxC,  max_regla_cat)
        
      } else { # DNF RULES (FALTA EL TRATAMIENTO DE VARIABLES CATEGORICAS)
        
        
        
        
        valNum <- mapply(FUN = ':', (max_regla_num + 1), (max_regla_num + nLabels), SIMPLIFY = FALSE)  
        regla_num <- lapply(X = valNum, FUN = function(x, rule) rule[x], regla )
        
        if(length(num_particip) > 0){
          
          fuzzy_sets <- dataset[["fuzzySets"]][1:nLabels, 1:3, num_particip, drop = F]
          crispSets <- dataset[["crispSets"]][1:nLabels, 1:2, num_particip, drop = F]
          
          #  Gets values for xmin, xmedio and xmax for fuzzy computation. 
          # The format is a matrix, which columns has at first value the number of numerical 
          # variable, and then, the values for xmin, xmedio, xmax, and only for values that participate in the rule
          n_matrices <- dim(fuzzy_sets)[3] 
          valuesFuzzy <- .getFuzzyValues(regla_num = regla_num, fuzzy = fuzzy_sets)
         
          #Gets values for xmin and xmax for crisp computation
          n_matricesCrisp <- dim(crispSets)[3]  
          valuesCrisp <- .getFuzzyValues(regla_num = regla_num, fuzzy = crispSets, crisp = TRUE)
        }
        
        #gr_perts <- lapply(X = noClass, FUN = .comparaDNF3, regla = regla, regla_num, cat_particip, num_particip,  max_regla_cat, max_regla_num, nLabels, fuzzy_sets, crispSets, valuesFuzzy, valuesCrisp)
        gr_perts <- .comparaDNF4(ejemplo = noClass, regla = regla, regla_num, cat_particip, num_particip,  max_regla_cat, max_regla_num, nLabels, fuzzy_sets, crispSets, valuesFuzzy, valuesCrisp)
        
        
        #gr_perts <- unlist(gr_perts)
      }
      
      
      values <- .get_values6(gr_perts = gr_perts, nombre_clases = dataset[["class_names"]], dataset = dataset[["data"]], targetClass = targetClass, examples_perClass = dataset[["examplesPerClass"]],cov = dataset[["covered"]], Ns = dataset[["Ns"]], N_vars = n_Vars + 1, por_cubrir = por_cubrir, marcar = marcar, test = test, difuso = difuso)
      
      #Compute fitness
      if(! marcar){
      
        fitness <- 0
        if(is.function(Objetivos[[1]]) && Pesos[1] > 0){ 
          fitness <- fitness + (Objetivos[[1]](values) * Pesos[1])
        }
        if(is.function(Objetivos[[2]]) && Pesos[2] > 0){ 
          fitness <- fitness + (Objetivos[[2]](values) * Pesos[2])
        }
        if(is.function(Objetivos[[3]]) && Pesos[3] > 0) {
          fitness <- fitness + (Objetivos[[3]](values) * Pesos[3])
        }      
        fitness <- fitness / (sum(Pesos))
        # cat("Ns:", values[[4]], " - Local Support: ", .LocalSupport(values), " - .confianza:", .confianza(values), " - Support: ", .Csupport(values)," - .coverage:", .coverage(values), " - Fitness: ", fitness, file = "", fill = TRUE)
        
        fitness #Return
      } else {
        
        values #Return
      }
      
    } else{
      0 #Return
    }
    
  } else {
    0 #Return
  }
  
}


.fitnessMESDIF <- function(regla, dataset, noClass, targetClass, por_cubrir, n_Vars, nLabels, max_regla, marcar = FALSE, Objetivos = c(.LocalSupport, .confianza, NULL, FALSE), Pesos = c(0.7,0.3,0), DNFRules = FALSE, difuso = FALSE ,test = FALSE, cate, num, NMEEF = FALSE){
  
  if( ! any(is.na(regla))) { #Si la regla no tiene NA se puede evaluar
    
    regla <- as.numeric(regla)
    participantes <- logical(length(max_regla))
    participantes <- .getParticipantes(regla = regla, max_regla = max_regla, DNFRules = DNFRules)
    
    
    #If it's not the empty rule
    if(any(participantes)){
      
      cat_particip <- which(cate & participantes)
      num_particip <- which(num & participantes)
      
      max_regla_cat <- max_regla[cat_particip]
      max_regla_num <- max_regla[num_particip]
      
      if(!DNFRules) { # CAN RULES
        
        #Split into numerical variables and categorical ones. (And participate in the rule)
        
        rule_cat <- regla[cat_particip]
        rule_num <- regla[num_particip]
        
        if(length(num_particip) > 0){
          fuzzy_sets <- dataset[["fuzzySets"]][1:nLabels, 1:3, num_particip, drop = F]
          crispSets <- dataset[["crispSets"]][1:nLabels, 1:2, num_particip, drop = F]
          #  Get values for xmin, xmedio and xmax for fuzzy computation.   
          n_matrices <- dim(fuzzy_sets)[3]  
          xmin <- fuzzy_sets[cbind(rule_num + 1, 1, seq_len(n_matrices))]
          xmax <- fuzzy_sets[cbind(rule_num + 1, 3, seq_len(n_matrices))]
          xmedio <- fuzzy_sets[cbind(rule_num + 1, 2, seq_len(n_matrices))]
          
          #Get values for xmin and xmax for crisp computation
          n_matricesCrisp <- dim(crispSets)[3]  
          xminC <- crispSets[cbind(rule_num + 1, 1, seq_len(n_matricesCrisp))]
          xmaxC <- crispSets[cbind(rule_num + 1, 2, seq_len(n_matricesCrisp))]
        }
        
        gr_perts <- .compara_CAN9(ejemplo = noClass, rule_cat = rule_cat, rule_num = rule_num, catParticip = cat_particip, numParticip = num_particip, xmin = xmin, xmedio = xmedio, xmax = xmax, n_matrices = n_matrices, xminCrisp = xminC, xmaxCrisp = xmaxC,  max_regla_cat)
        
      } else { # DNF RULES
        
        
        
        
        valNum <- mapply(FUN = ':', (max_regla_num + 1), (max_regla_num + nLabels), SIMPLIFY = FALSE)  
        regla_num <- lapply(X = valNum, FUN = function(x, rule) rule[x], regla )
        
        if(length(num_particip) > 0){
          
          fuzzy_sets <- dataset[["fuzzySets"]][1:nLabels, 1:3, num_particip, drop = F]
          crispSets <- dataset[["crispSets"]][1:nLabels, 1:2, num_particip, drop = F]
          
          #  Gets values for xmin, xmedio and xmax for fuzzy computation. 
          # The format is a matrix, which columns has at first value the number of numerical 
          # variable, and then, the values for xmin, xmedio, xmax, and only for values that participate in the rule
          n_matrices <- dim(fuzzy_sets)[3] 
          valuesFuzzy <- .getFuzzyValues(regla_num = regla_num, fuzzy = fuzzy_sets)
          
          #Gets values for xmin and xmax for crisp computation
          n_matricesCrisp <- dim(crispSets)[3]  
          valuesCrisp <- .getFuzzyValues(regla_num = regla_num, fuzzy = crispSets, crisp = TRUE)
        }
        
        gr_perts <- .comparaDNF4(ejemplo = noClass, regla = regla, regla_num, cat_particip, num_particip,  max_regla_cat, max_regla_num, nLabels, fuzzy_sets, crispSets, valuesFuzzy, valuesCrisp)
        

      }
      
      
      if(!DNFRules)
        values <- .get_values6(gr_perts = gr_perts, nombre_clases = dataset[["class_names"]], dataset = dataset[["data"]], targetClass = targetClass, examples_perClass = dataset[["examplesPerClass"]],cov = dataset[["covered"]], Ns = dataset[["Ns"]], N_vars = n_Vars + 1, por_cubrir = por_cubrir, marcar = marcar, test = test, difuso = difuso, NMEEF)
      else
        values <- .get_values6(gr_perts = gr_perts, nombre_clases = dataset[["class_names"]], dataset = dataset[["data"]], targetClass = targetClass, examples_perClass = dataset[["examplesPerClass"]],cov = dataset[["covered"]], Ns = dataset[["Ns"]], N_vars = n_Vars + 1, por_cubrir = por_cubrir, marcar = marcar, test = test, difuso = difuso, NMEEF)
      
      #Compute fitness
      if(! marcar){
        fitness <- numeric(4)
        fitness[1] <- if(is.function( Objetivos[[1]])) Objetivos[[1]](values) else 0
        fitness[2] <- if(is.function( Objetivos[[2]])) Objetivos[[2]](values) else 0
        fitness[3] <- if(is.function( Objetivos[[3]])) Objetivos[[3]](values) else 0
        if(! NMEEF)
          fitness #Return
        else 
          list(fit = fitness, covered = values[[13]]) # Return
      } else {
        
        values #Return
      }
      
    } else{
      c(0,0,0,0) #Return
    }
    
  } else {
    c(0,0,0,0) #Return
  }
  
}


#'
#' Obtains the belonging degree of every example of a dataset to a given rule
#' 
#' @param regla The rule to compare example. This rule must be in canonica vector representation. (See Rule.toRuleCANRepresentation function)
#' @param dataset The complete keel dataset object to get the examples
#' @param noClass a matrix with all examples without the class attribute. One examples PER COLUMN
#' @param nLabels number of fuzzy Labels that have numerical attributes
#' @param max_regla maximum value of all attributes ($conjuntos of the keel dataset)
#' @param cate logical vector indicating which attributes are categorical
#' @param num logical vector indicating which attributes are numerical
#' @param The T-norm to use. 0 to Minimum T-norm, 1 to Product T-norm.
#'
#' @return a numeric vector with the belonging degree of every example to the given rule.
#' 
.fitnessFuGePSD <- function(regla, dataset, noClass, nLabels, max_regla, cate, num, t_norm){
  
  
  if( ! any(is.na(regla))) { #Si la regla no tiene NA se puede evaluar
    
    regla <- as.integer(regla)
    participantes <- logical(length(max_regla))
    participantes <- .getParticipantes(regla = regla, max_regla = max_regla, DNFRules = FALSE)
    
    
    #If it's not the empty rule
    if(any(participantes)){
      
      cat_particip <- which(cate & participantes)
      num_particip <- which(num & participantes)
      
      max_regla_cat <- max_regla[cat_particip]
      max_regla_num <- max_regla[num_particip]
      
     
        
        #Split into numerical variables and categorical ones. (And participate in the rule)
        if(length(cat_particip) > 0){
          rule_cat <- regla[cat_particip]
        }
        
        if(length(num_particip) > 0){
          rule_num <- regla[num_particip]
          
          fuzzy_sets <- dataset[["fuzzySets"]][1:nLabels, 1:3, num_particip, drop = F]
          crispSets <- dataset[["crispSets"]][1:nLabels, 1:2, num_particip, drop = F]
          #  Get values for xmin, xmedio and xmax for fuzzy computation.   
          n_matrices <- dim(fuzzy_sets)[3]  
          
          xmin <- fuzzy_sets[cbind(rule_num + 1, 1, seq_len(n_matrices))]
          xmax <- fuzzy_sets[cbind(rule_num + 1, 3, seq_len(n_matrices))]
          xmedio <- fuzzy_sets[cbind(rule_num + 1, 2, seq_len(n_matrices))]
          
          #Get values for xmin and xmax for crisp computation
          n_matricesCrisp <- dim(crispSets)[3]  
          xminC <- crispSets[cbind(rule_num + 1, 1, seq_len(n_matricesCrisp))]
          xmaxC <- crispSets[cbind(rule_num + 1, 2, seq_len(n_matricesCrisp))]
        }
        
      #return
        Rule.compatibility(ejemplo = noClass, rule_cat = rule_cat, rule_num = rule_num, catParticip = cat_particip, numParticip = num_particip, xmin = xmin, xmedio = xmedio, xmax = xmax, n_matrices = n_matrices, max_cat = max_regla_cat, max_num = max_regla_num, t_norm = t_norm)
        
      
    }
  }
}
#         values <- .get_values6(gr_perts = gr_perts, nombre_clases = dataset[["class_names"]], dataset = dataset[["data"]], targetClass = targetClass, examples_perClass = dataset[["examplesPerClass"]],cov = dataset[["covered"]], Ns = dataset[["Ns"]], N_vars = n_Vars + 1, por_cubrir = por_cubrir, marcar = marcar, test = test, difuso = difuso)
#       
#       #Compute fitness
#       if(! marcar){
#         
#         fitness <- 0
#         if(is.function(Objetivos[[1]]) && Pesos[1] > 0){ 
#           fitness <- fitness + (Objetivos[[1]](values) * Pesos[1])
#         }
#         if(is.function(Objetivos[[2]]) && Pesos[2] > 0){ 
#           fitness <- fitness + (Objetivos[[2]](values) * Pesos[2])
#         }
#         if(is.function(Objetivos[[3]]) && Pesos[3] > 0) {
#           fitness <- fitness + (Objetivos[[3]](values) * Pesos[3])
#         }      
#         fitness <- fitness / (sum(Pesos))
#         # cat("Ns:", values[[4]], " - Local Support: ", .LocalSupport(values), " - .confianza:", .confianza(values), " - Support: ", .Csupport(values)," - .coverage:", .coverage(values), " - Fitness: ", fitness, file = "", fill = TRUE)
#         
#         fitness #Return
#       } else {
#         
#         values #Return
#       }
#       
#     } else{
#       0 #Return
#     }
#     
#   } else {
#     0 #Return
#   }
  







