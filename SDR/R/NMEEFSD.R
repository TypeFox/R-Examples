
#
#
#    Re-initialization operator based on coverage.
#  
#
#

.reInitPob <- function(elitePop, fitnessElite, coveredElite, crowdingDistance, pctVariables, cubiertoActual, dataset, maxRegla, cate, num, crispSets, targetClass, popSize){
  salida <- matrix(ncol = NCOL(elitePop), nrow = popSize)
  fitnessSalida <- matrix(ncol = 4, nrow = popSize*2)
  crowding <- numeric(popSize * 2)
  numVariables <- round(NCOL(elitePop) * pctVariables)
  coveredSalida <- matrix(FALSE, nrow = NROW(coveredElite), ncol = popSize)
  
  #Remove duplicated in elitePop and add this into salida
  noDuplicados <- which(! duplicated(elitePop))

  elitePop <- elitePop[noDuplicados,, drop = F]
  fitnessElite <- fitnessElite[noDuplicados,, drop = F]
  crowdingDistance <- crowdingDistance[noDuplicados]
  coveredElite <- coveredElite[, noDuplicados, drop = F]
  numIndividuos <- NROW(elitePop)

    salida[seq_len(numIndividuos), ] <- elitePop
    fitnessSalida[seq_len(NROW(fitnessElite)), ] <- fitnessElite
    crowding[seq_len(length(crowdingDistance))] <- crowdingDistance
    coveredSalida[,seq_len(NCOL(coveredElite))] <- coveredElite

  restantes <- popSize - numIndividuos
  # Initialization based on coverage of the rest of the population 
  # Find an example not covered and has the target class

  regla <- numeric(length(maxRegla))
  
  for(i in seq_len(restantes)){
    regla[] <- 0
    regla <- regla + maxRegla
    
    # Find an example not covered and has the target class
    noCubiertos <- which(! cubiertoActual)
    noCubiertos <- which( dataset[NROW(dataset), noCubiertos] == targetClass )
    
    if(length(noCubiertos) > 0)
           if(length(noCubiertos > 1))
            ejemplo <- dataset[seq_len(NROW(dataset) - 1), sample(noCubiertos, size = 1)]
           else
             ejemplo <- dataset[seq_len(NROW(dataset) - 1), noCubiertos]
         else
           ejemplo <- dataset[seq_len(NROW(dataset) - 1), sample(NCOL(dataset), size = 1)]
      
         #Get the values of all variable that cover the example
            valores <- numeric(length(ejemplo))
            for(ii in seq_len(length(valores))){
              if(cate[ii])
                valores[ii] <- ejemplo[ii] + 1
              
              if(num[ii]){ # Numerical variable
                #Get the interval the variable belongs to
               valores[ii] <- which( ejemplo[ii] > (crispSets[,1,ii] + .tolerance) & ejemplo[ii] <= (crispSets[,2,ii] + .tolerance)) - 1
              }
            }
            
            #Get the value of the selected variables that covers the example
            numInterv <- sample(numVariables, size = 1)
            vars <- sample(length(maxRegla), size = numInterv, replace = FALSE)
            
            regla[vars] <- valores[vars]
            numIndividuos <- numIndividuos + 1
            salida[numIndividuos,] <- regla
  }
  
  list(pop = salida, fitness = fitnessSalida, crowd = crowding, cov = coveredSalida)
}





#
#
# Selection operator for NMEEF-SD
#
#

.selectionNMEEF <- function(pop, popSize, rank, crowding, fitness, covered) {
  salida <- matrix(nrow = popSize, ncol = ncol(pop))
  fitnessSalida <- matrix(nrow = popSize, ncol = 4)
  coveredSalida <- matrix(nrow = NROW(covered), ncol = popSize)
  mating <- matrix(NA, nrow = 2, ncol = popSize)
  
  equals <- seq_len(popSize) # Tournaments among the same individual isnt allowed
 while(length(equals > 0)){
    mating[,equals] <- matrix(sample(seq_len(popSize), size = length(equals) * 2, replace = TRUE) , nrow = 2)
    equals <- which(mating[1,] == mating[2,])
  }
  
  # First, we compare the individuals by his rank value (less is better)
  winners1 <- which(rank[mating[1,]] < rank[mating[2,]])
  winners2 <- which(rank[mating[2,]] < rank[mating[1,]])
  
  pos <- 1
  if (length(winners1) > 0) {
    salida[seq_len(length(winners1)),] <- pop[mating[1,winners1],]
    fitnessSalida[seq_len(length(winners1)),] <- fitness[mating[1,winners1],]
    coveredSalida[,seq_len(length(winners1))] <- covered[,mating[1,winners1]]
    pos <- pos + length(winners1)
  }
  if (length(winners2) > 0) {
    salida[pos:(pos+length(winners2) -1),] <- pop[mating[2,winners2],]
    fitnessSalida[pos:(pos+length(winners2) -1 ),] <- fitness[mating[2,winners2],]
    coveredSalida[,pos:(pos+length(winners2) -1 )] <- covered[,mating[2,winners2]]
    pos <- pos + length(winners2)
  }
  
  #If there are ties, we solve it by his crowding distance (more is better)
  iguales <- which(rank[mating[2,]] == rank[mating[1,]])
  if (length(iguales) > 0) {
    mating <- mating[,iguales, drop = F]
    
    winners1 <- which(crowding[mating[1,]] >= crowding[mating[2,]])
    winners2 <- which(crowding[mating[2,]] > crowding[mating[1,]])
    if (length(winners1) > 0) {
      salida[pos:(pos+length(winners1) - 1),] <- pop[mating[1,winners1],]
      fitnessSalida[pos:(pos+length(winners1) -1),] <- fitness[mating[1,winners1],]
      coveredSalida[,pos:(pos+length(winners1) -1 )] <- covered[,mating[1,winners1]]
      pos <- pos + length(winners1)
    }
    if (length(winners2) > 0) {
      salida[pos:(pos+length(winners2) -1),] <- pop[mating[2,winners2],]
      fitnessSalida[pos:(pos+length(winners2) -1),] <- fitness[mating[2,winners2],]
      coveredSalida[,pos:(pos+length(winners2) -1 )] <- covered[,mating[2,winners2]]
    }
  }
  
  #Returns the selected population, the fitness of this individuals and his covered examples. (crowding distance will be computed after)
  list(population = salida, fitness = fitnessSalida, cov = coveredSalida)
  
}







#
#
#  Function to calculate dominance between individuals.
#
#

.calculateDominance <- function(q, p, strictDominance){
  #We calculate if p domain q only.
  
  dominaP <- F
  dominaQ <- F
  
  if(strictDominance){
    dominaP <- any(p > q)
    dominaQ <- any(p < q)
  } else {
    dominaQ <- any(p < q)
    dominaP <- any(p >= q)
  }
  
  #return
  if(dominaQ == dominaP)
    return(0L)
  if(dominaQ)
    return(1L)
  if(dominaP)
    return(-1L)
}






#
#
#  This function calculate the crowding distance of an individual in a front
#
#

.crowdingDistance <- function(pop){
  if(is.vector(pop)) 
    size <- 1
  else
    size <- NROW(pop)
  
  distance <- numeric(size)
  measures <- numeric(size)
  num_measures <- seq_len(NCOL(pop))
  
  
  if(size <= 2){
    distance <- Inf
  } else {
    #For every measure
    for( i in num_measures){
      measures <- pop[,i]
      #Order measures
      indices <- order(measures, decreasing = FALSE)
      
      #/**/ REMOVE THIS AND USE ORDER
      #indices <- .qsort(measures, 1, length(measures), seq_len(length(measures)))
      #indices <- indices$indices
      #/**/
      
      #Set boundary individuals distance as infinity
      distance[c(indices[1], indices[size])] <- Inf
      
      #Compute distance of the rest.
      for(j in 2:(size - 1)){
        distancia <- measures[indices[j + 1]] - measures[indices[j - 1]]
        if(distancia != 0){
          distance[indices[j]] <- distance[indices[j]] + (distancia / (measures[indices[size]] - measures[indices[1]]) )
        }
      }
    }
  }
  #Return
  distance
}






#
#
# This function fills the population that participate in the next generation of the evolutionary process
#
#

.fillPopulation <- function(fronts, numFronts, fitness, coveredFrentes, popSize, nObjs){
  suma <- 0
  if(is.vector(fronts[[1]]))
    cols <- length(fronts[[1]])
  else
    cols <- NCOL(fronts[[1]])
  
  newPop <- matrix(nrow = popSize, ncol = cols)
  newRank <- numeric(popSize * 2) + Inf  #No ranked indivudals cant be selected
  distance <- numeric(popSize)
  fit <- matrix(ncol = 4, nrow = popSize)
  frente <- 1
  coveredbyInd <- matrix(FALSE, nrow = NROW(coveredFrentes[[1]]), ncol= popSize * 2)
  #Indicate if the last front introduced in the population fits perfectly or we have to make the ordering by crowding distance of the front
  FitPerfectly <- TRUE
  
  #for(i in seq_len(length(fronts))){
  for(i in seq_len(numFronts)){
    #If front fits in newPop, we introduce completely in it
    frente <- i
    if(! is.vector(fronts[[i]]))
      rows <- NROW(fronts[[i]])
    else 
      rows <- 1
    
    if( rows + suma <= popSize ){
      
      #Calculate crowding Distance
      distance[(suma+1):(rows + suma)] <- .crowdingDistance(fitness[[i]][, seq_len(nObjs)])
      #Add front in the new population, and update the rest of parameters.
      newPop[(suma+1):(rows + suma), ] <- fronts[[i]]
      newRank[(suma+1):(rows + suma)] <- i - 1
      fit[(suma+1):(rows + suma), ] <- fitness[[i]]
      coveredbyInd[, (suma+1):(rows + suma)] <- coveredFrentes[[i]]
      suma <- suma + rows
    } else {
      break
    }
  }
  
    #If suma is less than popSize, front "frente" doesnt fit completely, so we must order by crowding distance
    if(suma < popSize){
      FitPerfectly <- FALSE
      porRellenar <- popSize - suma
      #Calculate crowding distance
      distancia <- .crowdingDistance(fitness[[frente]])
      
      #Order by crowding distance and get the best individuals till the population is filled
      #orden <- order(distancia, decreasing = TRUE)[seq_len(porRellenar)]
      orden <- .qsort(distancia, 1, length(distancia), index = seq_len(length(distancia)))
      orden <- orden$indices[length(orden$vector):(length(orden$vector) - porRellenar + 1) ]
      
      #Fill the population and update parameters
      newPop[(suma + 1):popSize, ] <- fronts[[frente]][orden, , drop = F]
      distance[(suma + 1):popSize] <- distancia[orden]
      newRank[(suma+1):popSize] <- frente - 1
      fit[(suma + 1):popSize, ] <- fitness[[frente]][orden, , drop = F]
      coveredbyInd[, (suma + 1):popSize] <- coveredFrentes[[frente]][,orden,drop = F]
    }
  
  list(population = newPop, distancia = distance, fitness = fit, rank = newRank, cov = coveredbyInd, fits = FitPerfectly)
  
  
}






#
#
# Mutation operator for NMEEF-SD
#
#

.mutateNMEEF <- function(cromosoma, variable, max_valor_variables, DNF_Rule){
  
  # Estos pesos hay que quitarlos una vez SDIGA funcione correctamente.
  mutation_type <- sample(x = 1:2, size = 1, prob = c(6/11, 5/11))   # Tipo 1 - tipo 2 aleatoriamente
  
  
  if(! DNF_Rule){  #Reglas can-nicas
    if(mutation_type == 1L){
      
      cromosoma[variable] <- max_valor_variables[variable] #Se pone el valor de no participacion
      
    } else {  #Asigna valor aleatorio en la variable (NO incluye valor de eliminacion)
      
      value <- sample(x = 0:(max_valor_variables[variable] - 1), size = 1)
      cromosoma[variable] <- value 
      
    }
    
  } else { #Reglas DNF
    
    variable <- variable + 1
    rango <- (max_valor_variables[variable - 1] + 1):max_valor_variables[variable]
    
    
    if(mutation_type == 1){  #Valor de no participaci-n de la variable
      
      cromosoma[rango] <- 0
      
    } else {  #Asigna valor aleatorio en la variable
      
      cromosoma[rango] <- sample(x = 0:1 , size = length(rango), replace = TRUE)
      
    }
    
  }
  
  
  
  cromosoma  # Return
  
}












#' @title Non-dominated Multi-objective Evolutionary algorithm for Extracting Fuzzy rules in Subgroup Discovery (NMEEF-SD)
#' @description Perfoms a subgroup discovery task executing the algorithm NMEEF-SD
#'
#' @param paramFile The path of the parameters file. \code{NULL} If you want to use training and test \code{keel} variables
#' @param training A \code{keel} class variable with training data.
#' @param test A \code{keel} class variable with training data.
#' @param output character vector with the paths of where store information file, rules file and test quality measures file, respectively.
#' @param seed An integer to set the seed used for generate random numbers.
#' @param nLabels Number of fuzzy labels defined in the datasets.
#' @param nEval An integer for set the maximum number of evaluations in the evolutive process.
#' @param popLength An integer to set the number of individuals in the population.
#' @param crossProb Sets the crossover probability. A number in [0,1].
#' @param mutProb Sets the mutation probability. A number in [0,1].
#' @param RulesRep Representation used in the rules. "can" for canonical rules, "dnf" for DNF rules.
#' @param Obj1 Sets the Objective number 1. See \code{Objective values} for more information about the possible values.
#' @param Obj2 Sets the Objective number 2. See \code{Objective values} for more information about the possible values.
#' @param Obj3 Sets the Objective number 3. See \code{Objective values} for more information about the possible values.
#' @param minCnf Sets the minimum confidence that must have a rule in the Pareto front for being returned. A number in [0,1].
#' @param reInitCoverage Sets if the algorithm must perform the reinitialitation based on coverage when it is needed. A string with "yes" or "no".
#' @param porcCob Sets the maximum percentage of variables that participate in the rules genereted in the reinitialitation based on coverage. A number in [0,1]
#' @param StrictDominance Sets if the comparison between individuals must be done by strict dominance or not. A string with "yes" or "no".
#' @param targetVariable The name or index position of the target variable (or class). It must be a categorical one.
#' @param targetClass A string specifing the value the target variable. \code{null} for search for all possible values.
#' 
#' 
#' @details This function sets as target variable the last one that appear in the KEEL file. If you want 
#'     to change the target variable, you can use \link{changeTargetVariable} for this objective.  
#'     The target variable MUST be categorical, if it is not, throws an error.
#'     
#'     If you specify in \code{paramFile} something distintc to \code{NULL} the rest of the parameters are
#'     ignored and the algorithm tries to read the file specified. See "Parameters file structure" below 
#'     if you want to use a parameters file.
#' 
#' @section How does this algorithm work?:
#'     NMEEF-SD is a multiobjetctive genetic algorithm based on a NSGA-II approach. The algorithm
#'     first makes a selection based on binary tournament and save the individuals in a offspring population.
#'     Then, NMEEF-SD apply the genetic operators over individuals in offspring population
#'     
#'     For generate the population which participate in the next iteration of the evoluationary process
#'     NMEEF-SD calculate the dominance among all individuals (join main population and offspring) and then, apply the NSGA-II fast sort algorithm to order
#'     the population by fronts of dominance, the first front is the non-dominanted front (or Pareto), the second is 
#'     where the individuals dominated by one individual are, the thirt front dominated by two and so on.
#'     
#'     To promove diversity NMEEF-SD has a mechanism of reinitialization of the population based on coverage
#'     if the Pareto doesnt evolve during a 5%% of the total of evaluations.
#'     
#'     At the final of the evolutionary process, the algorithm returns only the individuals in the Pareto front
#'     which has a confidence greater than a minimum confidence level.
#'     
#' @section Parameters file structure:
#'   The \code{paramFile} argument points to a file which has the necesary parameters for NMEEF-SD works.
#'   This file \strong{must} be, at least, those parameters (separated by a carriage return):
#'   \itemize{
#'     \item \code{algorithm}  Specify the algorithm to execute. In this case. "NMEEFSD"
#'     \item \code{inputData}  Specify two paths of KEEL files for training and test. In case of specify only the name of the file, the path will be the working directory.
#'     \item \code{seed}  Sets the seed for the random number generator
#'     \item \code{nLabels}  Sets the number of fuzzy labels to create when reading the files
#'     \item \code{nEval}  Set the maximun number of \strong{evaluations of rules} for stop the genetic process
#'     \item \code{popLength}  Sets number of individuals of the main population
#'     \item \code{ReInitCob}  Sets if NMEEF-SD do the reinitialization based on coverage. Values: "yes" or "no"  
#'     \item \code{crossProb}  Crossover probability of the genetic algorithm. Value in [0,1]
#'     \item \code{mutProb}  Mutation probability of the genetic algorithm. Value in [0,1]
#'     \item \code{RulesRep}  Representation of each chromosome of the population. "can" for canonical representation. "dnf" for DNF representation.
#'     \item \code{porcCob}  Sets the maximum percentage of variables participe in a rule when doing the reinitialization based on coverage. Value in [0,1]
#'     \item \code{Obj1} Sets the objective number 1. 
#'     \item \code{Obj2} Sets the objective number 2. 
#'     \item \code{Obj3} Sets the objective number 3. 
#'     \item \code{minCnf} Minimum confidence for returning a rule of the Pareto. Value in [0,1] 
#'     \item \code{StrictDominance} Sets if the comparison of individuals when calculating dominance must be using strict dominance or not. Values: "yes" or "no"
#'     \item \code{targetClass}  Value of the target variable to search for subgroups. The target variable \strong{is always the last variable.}. Use \code{null} to search for every value of the target variable
#'   }
#'   
#'   An example of parameter file could be:
#'  \preformatted{
#' algorithm = NMEEFSD
#' inputData = "irisd-10-1tra.dat" "irisd-10-1tra.dat" "irisD-10-1tst.dat"
#' outputData = "irisD-10-1-INFO.txt" "irisD-10-1-Rules.txt" "irisD-10-1-TestMeasures.txt"
#' seed = 1
#' RulesRep = can
#' nLabels = 3
#' nEval = 500
#' popLength = 51
#' crossProb = 0.6
#' mutProb = 0.1
#' ReInitCob = yes
#' porcCob = 0.5
#' Obj1 = comp
#' Obj2 = unus
#' Obj3 = null
#' minCnf = 0.6
#' StrictDominance = yes
#' targetClass = Iris-setosa
#' }
#'   
#'   @section Objective values:
#'      You can use the following quality measures in the ObjX value of the parameter file using this values:
#'       \itemize{
#'         \item Unusualness -> \code{unus}
#'         \item Crisp Support -> \code{csup}
#'         \item Crisp Confidence -> \code{ccnf}
#'         \item Fuzzy Support -> \code{fsup}
#'         \item Fuzzy Confidence -> \code{fcnf}
#'         \item Coverage -> \code{cove}
#'         \item Significance -> \code{sign}
#'       }
#'     
#'     If you dont want to use a objetive value you must specify \code{null}
#'  
#'   
#' @return The algorithm shows in the console the following results:
#' \enumerate{
#'  \item The parameters used in the algorithm
#'  \item The rules generated.
#'  \item The quality measures for test of every rule and the global results.
#' 
#'     Also, the algorithms save those results in the files specified in the \code{output} parameter of the algorithm or 
#'     in the \code{outputData} parameter in the parameters file.
#' }
#' 
#' @references Carmona, C., Gonzalez, P., del Jesus, M., & Herrera, F. (2010). NMEEF-SD: Non-dominated Multi-objective Evolutionary algorithm for Extracting Fuzzy rules in Subgroup Discovery. 
#' @examples  
#'    NMEEF_SD(paramFile = NULL, 
#'                training = habermanTra, 
#'                test = habermanTst, 
#'                output = c("optionsFile.txt", "rulesFile.txt", "testQM.txt"),
#'                seed = 0, 
#'                nLabels = 3,
#'                nEval = 300, 
#'                popLength = 100, 
#'                mutProb = 0.1,
#'                crossProb = 0.6,
#'                RulesRep = "can",
#'                Obj1 = "CSUP",
#'                Obj2 = "CCNF",
#'                Obj3 = "null",
#'                minCnf = 0.6,
#'                reInitCoverage = "yes",
#'                porcCob = 0.5,
#'                StrictDominance = "yes",
#'                targetClass = "positive"
#'                )  
#' \dontrun{
#'       NMEEF_SD(paramFile = NULL, 
#'                training = habermanTra, 
#'                test = habermanTst, 
#'                output = c("optionsFile.txt", "rulesFile.txt", "testQM.txt"),
#'                seed = 0, 
#'                nLabels = 3,
#'                nEval = 300, 
#'                popLength = 100, 
#'                mutProb = 0.1,
#'                crossProb = 0.6,
#'                RulesRep = "can",
#'                Obj1 = "CSUP",
#'                Obj2 = "CCNF",
#'                Obj3 = "null",
#'                minCnf = 0.6,
#'                reInitCoverage = "yes",
#'                porcCob = 0.5,
#'                StrictDominance = "yes",
#'                targetClass = "null"
#'                )
#'      }
#'  @export       
NMEEF_SD <- function(paramFile = NULL, 
                     training = NULL, 
                     test = NULL, 
                     output = c("optionsFile.txt", "rulesFile.txt", "testQM.txt"),
                     seed = 0, 
                     nLabels = 3,
                     nEval = 10000, 
                     popLength = 100, 
                     mutProb = 0.1,
                     crossProb = 0.6,
                     RulesRep = "can",
                     Obj1 = "CSUP",
                     Obj2 = "CCNF",
                     Obj3 = "null",
                     minCnf = 0.6,
                     reInitCoverage = "yes",
                     porcCob = 0.5,
                     StrictDominance = "yes",
                     targetVariable = NA,
                     targetClass = "null"
                     )
{
 
  if(is.null(paramFile)){
    #Generate our "parameters file"
    if(class(training) != "keel" | class(test) != "keel")
      stop("'training' or 'test' parameters is not a KEEL class")
    
    if(is.null(training) | is.null(test)) 
      stop("Not provided a 'test' or 'training' file and neither a parameter file. Aborting...")
    
    if(training[[1]] != test[[1]] )
      stop("datasets ('training' and 'test') does not have the same relation name.")
    
    if(length(output) != 3 )
      stop("You must specify three files to save the results.")
    
    parametros <- list(seed = seed, 
                       algorithm = "NMEEFSD",
                       outputData = output,
                       nEval = nEval, 
                       popLength = popLength,
                       nLabels = nLabels,
                       mutProb = mutProb,
                       crossProb = crossProb,
                       RulesRep = RulesRep,
                       Obj1 = Obj1, 
                       Obj2 = Obj2,
                       Obj3 = Obj3,
                       minCnf = minCnf,
                       reInitPob = reInitCoverage,
                       porcCob = porcCob,
                       StrictDominance = StrictDominance,
                       targetClass = targetClass,
                       targetVariable = if(is.na(targetVariable)) training$atributeNames[length(training$atributeNames)] else targetVariable)
  } else {

    # Parametros --------------------------
    parametros <- .read.parametersFile2(file = paramFile)  # parametros del algoritmo
    
    if(parametros[[1]] != "NMEEFSD") stop("Parameters file has parameters for another algorithm, no for \"NMEEF-SD\"")
   
    training <- read.keel(file = parametros$inputData[1])   # training data 
    test <- read.keel(file = parametros$inputData[2])        # test data
   
  }
  if(is.na(parametros$targetVariable))
    parametros$targetVariable <- training$atributeNames[length(training$atributeNames)]
  #Change target variable if it is neccesary
  training <- changeTargetVariable(training, parametros$targetVariable)
  test <- changeTargetVariable(test, parametros$targetVariable)
  #Check if the last variable is categorical.
  if(training$atributeTypes[length(training$atributeTypes)] != 'c' | test$atributeTypes[length(test$atributeTypes)] != 'c')
    stop("Target variable is not categorical.")
  #Set the number of fuzzy labels
  training <- modifyFuzzyCrispIntervals(training, parametros$nLabels)
  training$conjuntos <- .dameConjuntos(data_types = training$atributeTypes, max = training$max, n_labels = parametros$nLabels)
  test <- modifyFuzzyCrispIntervals(test, parametros$nLabels)
  test$conjuntos <- .dameConjuntos(data_types = test$atributeTypes, max = test$max, n_labels = parametros$nLabels)
  #Set Covered
  #training$covered <- logical(training$Ns)
  test$covered <- logical(test$Ns)
  
  file.remove(parametros$outputData[which(file.exists(parametros$outputData))])
 
  if(tolower(parametros$RulesRep) == "can"){
    DNF = FALSE
  } else {
    DNF = TRUE
    vars <-  Reduce(f = '+', x = training[["conjuntos"]], accumulate = TRUE)
    vars <- vars[length(vars)]
  }
  
  .show_parameters(params = parametros, train = training, test = test)
  contador <- 0
  
  Objetivos <- .parseObjetives(parametros = parametros, "NMEEFSD", DNF)
  
  if(all(is.na(Objetivos[1:3]))) stop("No objective values selected. You must select, at least, one objective value. Aborting...")
 
   cate <- training[["atributeTypes"]][- length(training[["atributeTypes"]])] == 'c'
  num <- training[["atributeTypes"]][- length(training[["atributeTypes"]])] == 'r' | training[["atributeTypes"]][- length(training[["atributeTypes"]])] == 'e'
  
  
  #---------------------------------------------------
  
  
  #----- OBTENCION DE LAS REGLAS -------------------
  if(parametros$targetClass != "null"){ # Ejecuci-n para una clase
    cat("\n", "\n", "Searching rules for only one value of the target class...", "\n", "\n", file ="", fill = TRUE) 
    reglas <- .findRule(parametros$targetClass, "NMEEFSD", training, parametros, DNF, cate, num, Objetivos, as.numeric(parametros[["porcCob"]]), parametros[["StrictDominance"]] == "yes", parametros[["reInitPob"]] == "yes", minCnf)
    
    if(! is.null(unlist(reglas)) ){
      if(! DNF) 
        reglas <-  matrix(unlist(reglas), ncol =  training[["nVars"]] + 1 , byrow = TRUE)
      else 
        reglas <-  matrix(unlist(reglas), ncol = vars + 1 , byrow = TRUE)
      
      for(i in seq_len(NROW(reglas))){
        cat("GENERATED RULE", i,   file = "", sep = " ",fill = TRUE)
        cat("GENERATED RULE", i,   file = parametros$outputData[2], sep = " ",fill = TRUE, append = TRUE)
        .print.rule(rule = as.numeric( reglas[i, - NCOL(reglas)] ), max = training$conjuntos, names = training$atributeNames, consecuente = reglas[i, NCOL(reglas)], types = training$atributeTypes,fuzzySets = training$fuzzySets, categoricalValues = training$categoricalValues, DNF, rulesFile = parametros$outputData[2])
        cat("\n","\n",  file = "", sep = "",fill = TRUE)
        cat("\n",  file = parametros$outputData[2], sep = "",fill = TRUE, append = TRUE)
      }
    } else {
      cat("No rules found with a confidence greater than the specified", file = "", fill = TRUE)
      cat("No rules found with a confidence greater than the specified\n", file = parametros$outputData[2], append = TRUE)
      reglas <- numeric()
    }
    
  } else {  #Ejecucion para todas las clases
    
    cat("\n", "\n", "Searching rules for all values of the target class...", "\n", "\n", file ="", fill = TRUE)  
    
    if(Sys.info()[1] == "Windows")
      reglas <- lapply(X = training$class_names, FUN = .findRule, "NMEEFSD", training, parametros, DNF, cate, num, Objetivos, as.numeric(parametros[["porcCob"]]), parametros[["StrictDominance"]] == "yes", parametros[["reInitPob"]] == "yes", minCnf  )
    else
      reglas <- parallel::mclapply(X = training$class_names, FUN = .findRule, "NMEEFSD", training, parametros, DNF, cate, num, Objetivos, as.numeric(parametros[["porcCob"]]), parametros[["StrictDominance"]] == "yes", parametros[["reInitPob"]] == "yes"   , mc.cores = parallel::detectCores() - 1, minCnf)
    
    if(! is.null(unlist(reglas)) ){
      if(! DNF) 
        reglas <-  matrix(unlist(reglas), ncol =  training[["nVars"]] + 1 , byrow = TRUE)
      else 
        reglas <-  matrix(unlist(reglas), ncol = vars + 1 , byrow = TRUE)
      
      #Print Rules (In windows, mclapply doesnt work, so the rules are printed by ".findRule" function)
    #if(Sys.info()[1] != "Windows")
      for(i in seq_len(NROW(reglas))){
        cat("GENERATED RULE", i,   file = "", sep = " ",fill = TRUE)
        cat("GENERATED RULE", i,   file = parametros$outputData[2], sep = " ",fill = TRUE, append = TRUE)
        .print.rule(rule = as.numeric( reglas[i, - NCOL(reglas)] ), max = training$conjuntos, names = training$atributeNames, consecuente = reglas[i, NCOL(reglas)], types = training$atributeTypes,fuzzySets = training$fuzzySets, categoricalValues = training$categoricalValues, DNF, rulesFile = parametros$outputData[2])
        cat("\n","\n",  file = "", sep = "",fill = TRUE)
        cat("\n",  file = parametros$outputData[2], sep = "",fill = TRUE, append = TRUE)
      }
    } else {
      cat("No rules found with a confidence greater than the specified", file = "", fill = TRUE)
      cat("No rules found with a confidence greater than the specified\n", file = parametros$outputData[2], append = TRUE)
      reglas <- numeric()
    }
    
  }
  
  #---------------------------------------------------
  
  cat("\n", "\n", "Testing rules...", "\n", "\n", file = "", sep = " ", fill = TRUE)
  
  #--------  Testeo de las reglas --------------------
  if(length(reglas) > 0){
  sumNvars <- 0
  sumCov <- 0
  sumFsup <- 0
  sumCsup <- 0
  sumCconf <- 0
  sumFconf <- 0
  sumUnus <- 0
  sumSign <- 0
  sumAccu <- 0
  
  n_reglas <- NROW(reglas)
  for(i in seq_len(n_reglas)){
    val <- .probeRule2(rule = reglas[i, - NCOL(reglas)], testSet = test, targetClass = reglas[i, NCOL(reglas)], numRule = i, parametros = parametros, Objetivos = Objetivos, Pesos = c(0.7,0.3, 0), cate = cate, num = num, DNF = DNF)
    test[["covered"]] <- val[["covered"]]
    sumNvars <- sumNvars + val[["nVars"]]
    sumCov <- sumCov + val[["coverage"]]
    sumFsup <- sumFsup + val[["fsupport"]]
    sumCconf <- sumCconf + val[["cconfidence"]]
    sumFconf <- sumFconf + val[["fconfidence"]]
    sumUnus <- sumUnus + val[["unusualness"]]
    sumSign <- sumSign + val[["significance"]]
    sumAccu <- sumAccu + val[["accuracy"]]
  }
  
  
  
  #Medidas de calidad globales
  cat("Global:", file ="", fill = TRUE)
  cat(paste("\t - N_rules:", NROW(reglas), sep = " "),
      paste("\t - N_vars:", round(sumNvars / n_reglas, 6), sep = " "),
      paste("\t - Coverage:", round(sumCov / n_reglas, 6), sep = " "),
      paste("\t - Significance:", round(sumSign / n_reglas, 6), sep = " "),
      paste("\t - Unusualness:", round(sumUnus / n_reglas, 6), sep = " "),
      paste("\t - Accuracy:", round(sumAccu / n_reglas, 6), sep = " "),
      paste("\t - CSupport:", round(sum(test[["covered"]] / test[["Ns"]]), 6), sep = " "),
      paste("\t - FSupport:", round(sumFsup / n_reglas, 6), sep = " "),
      paste("\t - FConfidence:", round(sumFconf / n_reglas, 6), sep = " "),
      paste("\t - CConfidence:", round(sumCconf / n_reglas, 6), sep = " "),
      file = "", sep = "\n"
  )
  
  #Medidas de calidad globales (Save in testMeasures File FOR THE WEB INTERFACE)

  
  #Medidas de calidad globales (Save in testMeasures File)
  
  cat( "Global:",
       paste("\t - N_rules:", nrow(reglas), sep = " "),
       paste("\t - N_vars:", round(sumNvars / n_reglas, 6), sep = " "),
       paste("\t - Coverage:", round(sumCov / n_reglas, 6), sep = " "),
       paste("\t - Significance:", round(sumSign / n_reglas, 6), sep = " "),
       paste("\t - Unusualness:", round(sumUnus / n_reglas, 6), sep = " "),
       paste("\t - Accuracy:", round(sumAccu / n_reglas, 6), sep = " "),
       paste("\t - CSupport:", round(sum(test[["covered"]] / test[["Ns"]]), 6), sep = " "),
       paste("\t - FSupport:", round(sumFsup / n_reglas, 6), sep = " "),
       paste("\t - FConfidence:", round(sumFconf / n_reglas, 6), sep = " "),
       paste("\t - CConfidence:", round(sumCconf / n_reglas, 6), sep = " "),
       file = parametros$outputData[3], sep = "\n", append = TRUE
  )
  } else {
    cat("No rules for testing", file = "", fill = TRUE)
    cat("No rules for testing", file = parametros$outputData[2], append = TRUE)
  }
  #---------------------------------------------------
  
}