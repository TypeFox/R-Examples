




#-------------------------------------------------------------------------------------
#                Operador de mutaci-n
#
#  - cromosoma:                 cromosoma a modificar
#  - max_valor_variables:       vector con los valores m-ximos de las variables (aqu- el m-ximo indicar- valor de no participaci-n.)
#  - DNF_Rule:                  -Estoy usando reglas DNF-
#-------------------------------------------------------------------------------------

.mutate <- function(cromosoma, variable, max_valor_variables, DNF_Rule){
                                            
 
  mutation_type <- sample(x = 1:2, size = 1)   # Tipo 1 - tipo 2 aleatoriamente
  
  
  if(! DNF_Rule){  #Reglas can-nicas
    if(mutation_type == 1L){
      
      cromosoma[variable] <- max_valor_variables[variable] #Se pone el valor de no participacion
      
    } else {  #Asigna valor aleatorio en la variable (Incluye valor de eliminacion)
      
      value <- sample(x = 0:(max_valor_variables[variable] ), size = 1)
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











#' @title Subgroup Discovery Iterative Genetic Algorithm (SDIGA)
#' @description Perfoms a subgroup discovery task executing the algorithm SDIGA
#' 
#' @param parameters_file The path of the parameters file. \code{NULL} If you want to use training and test \code{keel} variables
#' @param training A \code{keel} class variable with training data.
#' @param test A \code{keel} class variable with training data.
#' @param output character vector with the paths of where store information file, rules file and test quality measures file, respectively.
#' @param seed An integer to set the seed used for generate random numbers.
#' @param nLabels Number of fuzzy labels defined in the datasets.
#' @param nEval An integer for set the maximum number of evaluations in the evolutive process.
#' @param popLength An integer to set the number of individuals in the population.
#' @param mutProb Sets the mutation probability. A number in [0,1].
#' @param RulesRep Representation used in the rules. "can" for canonical rules, "dnf" for DNF rules.
#' @param Obj1 Sets the Objective number 1. See \code{Objective values} for more information about the possible values.
#' @param w1 Sets the weight of \code{Obj1}.
#' @param Obj2 Sets the Objective number 2. See \code{Objective values} for more information about the possible values.
#' @param w2 Sets the weight of \code{Obj2}.
#' @param Obj3 Sets the Objective number 3. See \code{Objective values} for more information about the possible values.
#' @param w3 Sets the weight of \code{Obj3}.
#' @param minConf Sets the minimum confidence that must have the rule returned by the genetic algorithm after the local optimitation phase. A number in [0,1].
#' @param lSearch Sets if the local optimitation phase must be performed. A string with "yes" or "no".
#' @param targetVariable A string with the name or an integer with the index position of the target variable (or class). It must be a categorical one.
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
#'     This algorithm has a genetic algorithm in his core. This genetic algorithm returns only the best
#'     rule of the population and it is executed so many times until a stop condition is reached. The stop condition is 
#'     that the rule returned must cover at least one new example (not covered by previous rules) and must have a confidence
#'     greater than a minimum.
#'     
#'    After returning the rule, a local improvement could be applied for make the rule more general. This local improve is done
#'    by means of a hill-climbing local search.
#'    
#'    The genetic algorithm cross only the two best individuals. But the mutation operator is applied over all the population,
#'    individuals from cross too.
#'    
#' @section Parameters file structure:
#'   The \code{parameters_file} argument points to a file which has the necesary parameters for SDIGA works.
#'   This file \strong{must} be, at least, those parameters (separated by a carriage return):
#'   \itemize{
#'     \item \code{algorithm}  Specify the algorithm to execute. In this case. "SDIGA"
#'     \item \code{inputData}  Specify two paths of KEEL files for training and test. In case of specify only the name of the file, the path will be the working directory.
#'     \item \code{seed}  Sets the seed for the random number generator
#'     \item \code{nLabels}  Sets the number of fuzzy labels to create when reading the files
#'     \item \code{nEval}  Set the maximun number of \strong{evaluations of rules} for stop the genetic process
#'     \item \code{popLength}  Sets number of individuals of the main population
#'     \item \code{mutProb}  Mutation probability of the genetic algorithm. Value in [0,1]
#'     \item \code{RulesRep}  Representation of each chromosome of the population. "can" for canonical representation. "dnf" for DNF representation.
#'     \item \code{Obj1} Sets the objective number 1. 
#'     \item \code{w1} Sets the weigth assigned to the objective number 1. Value in [0,1]
#'     \item \code{Obj2} Sets the objective number 2. 
#'     \item \code{w2} Sets the weigth assigned to the objective number 2. Value in [0,1]
#'     \item \code{Obj3} Sets the objective number 3. 
#'     \item \code{w3} Sets the weigth assigned to the objective number 3. Value in [0,1]
#'     \item \code{minConf} Sets the minimum confidence of the rule for checking the stopping criteria of the iterative process
#'     \item \code{lSearch} Perform the local search algorithm after the execution of the genetic algorithm? Values: "yes" or "no"
#'     \item \code{targetClass}  Value of the target variable to search for subgroups. The target variable \strong{is always the last variable.}. Use \code{null} to search for every value of the target variable
#'   }
#'   
#'   An example of parameter file could be:
#'  \preformatted{
#' algorithm = SDIGA
#' inputData = "irisD-10-1tra.dat" "irisD-10-1tst.dat"
#' outputData = "irisD-10-1-INFO.txt" "irisD-10-1-Rules.txt" "irisD-10-1-TestMeasures.txt"
#' seed = 0
#' nLabels = 3
#' nEval = 500
#' popLength = 100
#' mutProb = 0.01
#' minConf = 0.6
#' RulesRep = can
#' Obj1 = Comp
#' Obj2 = Unus
#' Obj3 = null
#' w1 = 0.7
#' w2 = 0.3
#' w3 = 0.0
#' lSearch = yes
#' }
#' 
#' 
#' @section Objective values:
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
#' }
#' 
#'     Also, the algorithms save those results in the files specified in the \code{output} parameter of the algorithm or 
#'     in the \code{outputData} parameter in the parameters file.
#' 
#' 
#' 
#' @references 
#' M. J. del Jesus, P. Gonzalez, F. Herrera, and M. Mesonero, "Evolutionary
#' Fuzzy Rule Induction Process for Subgroup Discovery: A case study in
#' marketing," IEEE Transactions on Fuzzy Systems, vol. 15, no. 4, pp.
#' 578-592, 2007.
#' 
#' @examples 
#' SDIGA(parameters_file = NULL, 
#'       training = habermanTra, 
#'       test = habermanTst, 
#'       output = c("optionsFile.txt", "rulesFile.txt", "testQM.txt"),
#'       seed = 0, 
#'       nLabels = 3,
#'       nEval = 300, 
#'       popLength = 100, 
#'       mutProb = 0.01, 
#'       RulesRep = "can",
#'       Obj1 = "CSUP", 
#'       w1 = 0.7,
#'       Obj2 = "CCNF",
#'       w2 = 0.3,
#'       Obj3 = "null",
#'       w3 = 0,
#'       minConf = 0.6,
#'       lSearch = "yes",
#'       targetClass = "positive")
#' \dontrun{
#' SDIGA(parameters_file = NULL, 
#'       training = habermanTra, 
#'       test = habermanTst, 
#'       output = c("optionsFile.txt", "rulesFile.txt", "testQM.txt"),
#'       seed = 0, 
#'       nLabels = 3,
#'       nEval = 300, 
#'       popLength = 100, 
#'       mutProb = 0.01, 
#'       RulesRep = "can",
#'       Obj1 = "CSUP", 
#'       w1 = 0.7,
#'       Obj2 = "CCNF",
#'       w2 = 0.3,
#'       Obj3 = "null",
#'       w3 = 0,
#'       minConf = 0.6,
#'       lSearch = "yes",
#'       targetClass = "positive")
#'       }
#' @export
SDIGA <- function(parameters_file = NULL, 
                  training = NULL, 
                  test = NULL, 
                  output = c("optionsFile.txt", "rulesFile.txt", "testQM.txt"),
                  seed = 0, 
                  nLabels = 3,
                  nEval = 10000, 
                  popLength = 100, 
                  mutProb = 0.01, 
                  RulesRep = "can",
                  Obj1 = "CSUP", 
                  w1 = 0.7,
                  Obj2 = "CCNF",
                  w2 = 0.3,
                  Obj3 = "null",
                  w3 = 0,
                  minConf = 0.6,
                  lSearch = "yes",
                  targetVariable = NA,
                  targetClass = "null")
{
  
  
  
  if(is.null(parameters_file)){
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
                       algorithm = "SDIGA",
                       outputData = output,
                       nEval = nEval, 
                       popLength = popLength,
                       nLabels = nLabels,
                       mutProb = mutProb,
                       RulesRep = RulesRep,
                       Obj1 = Obj1, 
                       w1 = w1, 
                       Obj2 = Obj2, 
                       w2 = w2,
                       Obj3 = Obj3,
                       w3 = w3, 
                       lSearch = lSearch,
                       minConf = minConf,
                       targetClass = targetClass,
                       targetVariable = if(is.na(targetVariable)) training$atributeNames[length(training$atributeNames)] else targetVariable)
      
  } else {
  # Parametros --------------------------
  parametros <- .read.parametersFile2(file = parameters_file)  # parametros del algoritmo
    if(parametros$algorithm != "SDIGA")
      stop("Parameters file is not for \"SDIGA\"")
  
  test <- read.keel(file = parametros$inputData[2])        # test data
  
  training <- read.keel(file = parametros$inputData[1])   # training data
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
  training$covered <- logical(training$Ns)
  test$covered <- logical(test$Ns)
  
  file.remove(parametros$outputData[which(file.exists(parametros$outputData))])

  
  if(tolower(parametros$RulesRep) == "can"){
    DNF = FALSE
  } else {
    DNF = TRUE
  }
  
  Objetivos <- .parseObjetives(parametros = parametros, "SDIGA", DNF)
  if(all(is.na(Objetivos[1:3]))) stop("No objetive values selected. You must select, at least, one objective value. Aborting...")

  Pesos <- c(parametros$w1, parametros$w2, parametros$w3)
  if(sum(Pesos) == 0) stop("Sum of weigths must be a value greater than zero.")
  
  Mejor <- TRUE  
  
  reglas <- list()
  
  .show_parameters(params = parametros, train = training, test = test)
  contador <- 0
  
  cate <- training[["atributeTypes"]][- length(training[["atributeTypes"]])] == 'c'
  num <- training[["atributeTypes"]][- length(training[["atributeTypes"]])] == 'r' | training[["atributeTypes"]][- length(training[["atributeTypes"]])] == 'e'
  
  
  #---------------------------------------------------
  
  
  #----- OBTENCION DE LAS REGLAS -------------------
  if(parametros$targetClass != "null"){ # Ejecuci-n para una clase
    #Check if target class is valid
    targetClass <- parametros$targetClass
    if(! any(training$class_names == targetClass)) stop("No valid target value provided.")
    cat("\n", "\n", "Searching rules for only one value of the target class...", "\n", "\n", file ="", fill = TRUE)  
    
    
    cat(" - Target value:", targetClass , file = "", sep = " ", fill = TRUE)
    cat("\n - Target value:", targetClass , file = parametros$outputData[2], sep = " ", fill = TRUE, append = TRUE)
    primera_regla <- TRUE
    Mejor = TRUE
    por_cubrir = training$examplesPerClass[[targetClass]]
    
    while(Mejor){
      Mejor <- FALSE
      
      rule <- .ejecutarga(algorithm = "SDIGA", dataset = training, targetClass = targetClass, n_vars = training$nVars, por_cubrir = por_cubrir, nLabels = parametros$nLabels, N_evals = parametros$nEval,  tam_pob = parametros$popLength, p_mut = parametros$mutProb, seed = parametros$seed, Objetivos = Objetivos, Pesos = Pesos, DNFRules = DNF, cate = cate, num = num)
      maxRule <-  if(!DNF) training$conjuntos else c(0, Reduce(f = '+', x = training[["conjuntos"]], accumulate = TRUE))
      values <- .fit13(regla = rule, dataset = training, noClass = matrix(unlist(.separar(training)), nrow = length(training[[2]]) - 1, ncol = length(training[[7]])), targetClass = targetClass, por_cubrir = por_cubrir, n_Vars = training$nVars, nLabels = parametros$nLabels, max_regla = maxRule, marcar = TRUE, Objetivos = Objetivos, Pesos = Pesos, DNFRules = DNF, difuso = Objetivos[[4]], test = TRUE, cate = cate, num = num)[[2]]
      
      if(tolower(parametros$lSearch) == "yes") rule <- .Busqueda_Local(att_obj = targetClass, regla = rule, DNF_Rules = DNF, dataset = training , confianza_minima = parametros$minConf, x = values, max_regla = maxRule, por_cubrir = por_cubrir, nLabels = parametros$nLabels, Objetivos = Objetivos, cate = cate, num = num)
      x <- .marcar_ejemplos(regla = rule, dataset = training, targetClass = targetClass, nVars = training$nVars, maxRegla = maxRule, por_cubrir = por_cubrir, nLabels = parametros$nLabels, Objetivos = Objetivos, Pesos = Pesos, DNFRules = DNF, cate = cate, num = num)
      
      if(x$cubreNuevos && x$confidence > parametros$minConf || primera_regla){
        primera_regla <- FALSE
        Mejor <- TRUE
        contador <- contador + 1
        training$covered <- x$covered[[1]]
        cat("\n"," GENERATED RULE", contador, ":",file = "", sep = " ", fill = TRUE)
        cat("\n"," GENERATED RULE", contador, ":",file = parametros$outputData[2], sep = " ", fill = TRUE, append = TRUE)
        .print.rule(rule = rule, max = training$conjuntos, names = training$atributeNames, consecuente = targetClass, types = training$atributeTypes,fuzzySets = training$fuzzySets, categoricalValues = training$categoricalValues, DNF, rulesFile = parametros$outputData[2])
        cat("\n", file = "")
        rule[length(rule) + 1] <- targetClass
        reglas[[contador]] <- rule
        por_cubrir <- x$porCubrir
        if(por_cubrir <= 0) Mejor <- FALSE #Si no quedan ejemplos por cubrir no volvemos a llamar al gen-tico.
        
      } else {
        cat(" GENERATED RULE", ":",file = "", sep = " ", fill = TRUE)
        cat("# Invalid (Low confidence or support)", "\n","\n", file = "", sep= " ", fill = TRUE)
        
        cat("\n GENERATED RULE", ":", "\n",
          "# Invalid (Low confidence or support)", "\n", file = parametros$outputData[2], sep= " ", fill = TRUE, append = TRUE)
        
      }
      
    }
    
    
    
  } else {  #Ejecucion para todas las clases
    cat("\n", "\n", "Searching rules for all values of the target class...", "\n", "\n", file ="", fill = TRUE)  
    for(i in seq_len(length(training[["class_names"]]))) {
    #for(i in training$class_names[ seq_len(length(training$class_names)) ]){
      targetClass <- training[["class_names"]][i]
      cat(" - Target value:", targetClass , file = "", sep = " ", fill = TRUE)
      cat(" \n - Target value:", targetClass , file = parametros$outputData[2], sep = " ", fill = TRUE, append = TRUE)
      primera_regla <- TRUE
      Mejor = TRUE
      por_cubrir = training$examplesPerClass[[i]]
      
      while(Mejor){
        Mejor <- FALSE
        
        rule <- .ejecutarga(algorithm = "SDIGA", dataset = training, targetClass = targetClass, n_vars = training$nVars, por_cubrir = por_cubrir, nLabels = parametros$nLabels, N_evals = parametros$nEval,  tam_pob = parametros$popLength, p_mut = parametros$mutProb, seed = parametros$seed, Objetivos = Objetivos, Pesos = Pesos, DNFRules = DNF, cate = cate, num = num)
        maxRule <-  if(!DNF) training$conjuntos else c(0, Reduce(f = '+', x = training[["conjuntos"]], accumulate = TRUE))
      
        values <- .fit13(regla = rule, dataset = training, noClass = matrix(unlist(.separar(training)), nrow = length(training[[2]]) - 1, ncol = length(training[[7]])), targetClass = targetClass, por_cubrir = por_cubrir, n_Vars = training$nVars, nLabels = parametros$nLabels, max_regla = maxRule, marcar = TRUE, Objetivos = Objetivos, Pesos = Pesos, DNFRules = DNF, difuso = Objetivos[[4]], test = TRUE, cate = cate, num = num)[[2]]
        
        if(tolower(parametros$lSearch) == "yes"){
          rule <- .Busqueda_Local(att_obj = targetClass, regla = rule, DNF_Rules = DNF, dataset = training , confianza_minima = parametros$minConf, x = values, max_regla = maxRule, por_cubrir = por_cubrir, nLabels = parametros$nLabels, Objetivos = Objetivos, cate = cate, num = num)
        }
       
          x <- .marcar_ejemplos(regla = rule, dataset = training, targetClass = targetClass, nVars = training$nVars, maxRegla = maxRule, por_cubrir = por_cubrir, nLabels = parametros$nLabels, Objetivos = Objetivos, Pesos = Pesos, DNFRules = DNF, cate = cate, num = num)
        
        if(x$cubreNuevos && x$confidence > parametros$minConf || primera_regla){
          primera_regla <- FALSE
          Mejor <- TRUE
          contador <- contador + 1
          training$covered <- x$covered[[1]]
          cat("\n"," GENERATED RULE", contador, ":",file = "", sep = " ", fill = TRUE)
          cat("\n"," GENERATED RULE", contador, ":",file = parametros$outputData[2], sep = " ", fill = TRUE, append = TRUE)
          .print.rule(rule = rule, max = training$conjuntos, names = training$atributeNames, consecuente = targetClass, types = training$atributeTypes,fuzzySets = training$fuzzySets, categoricalValues = training$categoricalValues, DNF, rulesFile = parametros$outputData[2])
          cat("\n", file = "")
          rule[length(rule) + 1] <- targetClass
          reglas[[contador]] <- rule
          por_cubrir <- x$porCubrir
          if(por_cubrir <= 0) Mejor <- FALSE #Si no quedan ejemplos por cubrir no volvemos a llamar al gen-tico.
          
        } else {
          cat(" GENERATED RULE", ":",file = "", sep = " ", fill = TRUE)
          cat("# Invalid (Low confidence or support)", "\n","\n", file = "", sep= " ", fill = TRUE)
          
          cat("\n GENERATED RULE", ":", "\n",
              "# Invalid (Low confidence or support)", "\n","\n", file = parametros$outputData[2], sep= " ", fill = TRUE, append = TRUE)
          
        }
        
      }
    }
    
    
  }
  
  #---------------------------------------------------
  
  cat("\n", "\n", "Testing rules...", "\n", "\n", file = "", sep = " ", fill = TRUE)
  
  #--------  Testeo de las reglas --------------------
  sumNvars <- 0
  sumCov <- 0
  sumFsup <- 0
  sumCsup <- 0
  sumCconf <- 0
  sumFconf <- 0
  sumUnus <- 0
  sumSign <- 0
  sumAccu <- 0
  
  n_reglas <- length(reglas)
  for(i in 1:n_reglas){
    val <- .probeRule2(rule = reglas[[i]][-length(reglas[[i]])], testSet = test, targetClass = reglas[[i]][length(reglas[[i]])], numRule = i, parametros = parametros, Objetivos = Objetivos, Pesos = Pesos, cate = cate, num = num, DNF = DNF)
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
  cat(paste("\t - N_rules:", length(reglas), sep = " "),
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
  

  cat( "Global:",
     paste("\t - N_rules:", length(reglas), sep = " "),
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
  
  
  
  
  #---------------------------------------------------
  
  #reglas  # Return
  
}



.probeRule2 <- function(rule, testSet, targetClass, numRule, parametros, Objetivos, Pesos, cate, num, DNF = FALSE){
  stopifnot(class(testSet) == "keel")
   
  maxRegla <- .dameConjuntos(data_types = testSet[[3]], max = testSet[[5]], n_labels = parametros$nLabels)
   
  if(DNF) maxRegla <- c(0,Reduce(f= '+', x = maxRegla, accumulate = TRUE))
  
  # OJO QUE CUENTAS LAS REGLAS ANTERIOREs, LOS EJEMPLOS POR CUBRIR DE CADA REGLA NO SON LOS INICIALES
  # HAY QUE SUBIR LA VARIABLE POR_CUBRIR A UN NIVEL SUPERIOR
  p <- .fit13(regla = rule, dataset = testSet, noClass = matrix(unlist(.separar(testSet)), nrow = length(cate)), targetClass = targetClass, por_cubrir = testSet$examplesPerClass[[targetClass]], n_Vars = testSet$nVars, nLabels = parametros$nLabels, max_regla = maxRegla, marcar = TRUE, Objetivos = Objetivos, Pesos = Pesos, DNFRules = DNF, difuso = Objetivos[[4]], test = TRUE, cate = cate, num = num)
  values <- p[[2]]
  testSet[["covered"]] <- testSet[["covered"]] | p[[1]] #For the global quality measure
  
  Cov <- round(.coverage(values), 6)
  sig <- round(.significance(values), 6)
  unus <- round(.unusualness(values), 6)
  acc <- round(.accuracy(values), 6)
  Csup <- round(.Csupport(values), 6)
  Fsup <- round(.Fsupport(values), 6)
  Ccnf <- round(.confianza(values), 6)
  Fcnf <- round(.confianzaDifusa(values), 6)
  
  if(DNF) {
    participantes <- .getParticipantes(regla = rule,  max_regla = maxRegla, DNFRules = TRUE) 
    nVars <- sum(participantes) + 1
  } else{
    nVars <- sum(rule < testSet[["conjuntos"]]) + 1 # +1 Por la variable del consecuente.
  }
  cat("Rule", numRule,":", file = "", sep = " ", fill = TRUE)
  cat(paste("\t - N_vars:", nVars, sep = " "),
      paste("\t - Coverage:", Cov, sep = " "),
      paste("\t - Significance:", sig, sep = " "),
      paste("\t - Unusualness:", unus, sep = " "),
      paste("\t - Accuracy:", acc, sep = " "),
      paste("\t - CSupport:", Csup, sep = " "),
      paste("\t - FSupport:", Fsup, sep = " "),
      paste("\t - CConfidence:", Ccnf, sep = " "),
      paste("\t - FConfidence:", Fcnf, sep = " "),
      file = "", sep = "\n"
  )
  cat("\n", file = "", sep = "\n")
  
  
  
  
  cat(paste("Rule", numRule,":"),
      paste("\t - N_vars:", nVars, sep = " "),
      paste("\t - Coverage:", Cov, sep = " "),
      paste("\t - Significance:", sig, sep = " "),
      paste("\t - Unusualness:", unus, sep = " "),
      paste("\t - Accuracy:", acc, sep = " "),
      paste("\t - CSupport:", Csup, sep = " "),
      paste("\t - FSupport:", Fsup, sep = " "),
      paste("\t - CConfidence:", Ccnf, sep = " "),
      paste("\t - FConfidence:", Fcnf, "\n", sep = " "),
      file = parametros$outputData[3], sep = "\n", append = TRUE
  )
  
#Return
    list( covered = testSet[["covered"]], 
                nVars = nVars,
                coverage = Cov, 
                significance = sig, 
                unusualness = unus,
                accuracy = acc,
                csupport = Csup,
                fsupport = Fsup,
                cconfidence = Ccnf,
                fconfidence = Fcnf) 
}



#--------------------------------------------------------------------------------------------------
#                 Borra una variable de una regla
#
# - regla: La regla a modificar
# - variable: La variable a eliminar
# - max_valor_variables: Numero de conjuntos difusos para cada variable
# - DNF_Rules: -Uso reglas DNF-
#
#--------------------------------------------------------------------------------------------------

.borrar_gen <- function(regla, variable, max_valor_variables, DNF_Rules){
  
  if(!DNF_Rules){ # Reglas canonicas
    
    regla[variable] <- max_valor_variables[variable]  #valor de no participacion
    
    
  } else {   #Reglas DNF
    
   rango <- (max_valor_variables[variable] + 1):max_valor_variables[variable + 1]
   regla[rango] <- 0
   
  }
  
  regla   #Return
  
}




#--------------------------------------------------------------------------------------------------
#         B-squeda local etapa de post-procesamiento de SDIGA
#
# - regla: La regla a optimizar
# - DNF_Rules: -Uso reglas DNF-
# - dataset: el conjunto de ejemplos marcados en caso de est-n cubiertos.
# - max_valor_variables: n-mero de conjuntos difusos de cada variable
# - .confianza_minima: valor m--nimo de .confianza a dar
# - Valores devueltos por la funcion ejemplos_cubiertos
#
#--------------------------------------------------------------------------------------------------


.Busqueda_Local <- function(att_obj, regla, DNF_Rules, dataset, confianza_minima, x, max_regla, por_cubrir, nLabels, Objetivos, cate, num){
  mejor_regla <- regla
   soporte_regla <- .significance(x)
  participantes <- .getParticipantes(regla = regla, max_regla = max_regla, DNFRules = DNF_Rules)
  
  if(soporte_regla == 1 || sum(participantes) == 1 ){
    return(regla) # If local support it is 1 or rule has only one attribute, we can not improve the rule
  }
  
    mejor_soporte <-  soporte_regla
    mejor = TRUE
    longitud = if(DNF_Rules) length(max_regla) - 1 else length(max_regla)
    while(mejor){
      mejor = FALSE
      participantes <- .getParticipantes(regla = mejor_regla, max_regla = max_regla, DNFRules = DNF_Rules)
      for(i in seq_len(longitud) ){ # Para cada gen de la regla
     
        if(participantes[i]){
          regla_m <- .borrar_gen(regla = mejor_regla, variable = i, max_valor_variables = max_regla, DNF_Rules = DNF_Rules)
          x1 <- .fit13(regla = regla_m, dataset = dataset, noClass = matrix(unlist(.separar(dataset)), nrow = length(dataset[[2]]) - 1, ncol = length(dataset[[7]])), targetClass = att_obj, por_cubrir = por_cubrir, n_Vars = dataset$nVars, nLabels = nLabels, max_regla = max_regla, marcar = TRUE, Objetivos = Objetivos, DNFRules = DNF_Rules, difuso = Objetivos[[4]], test = TRUE, cate = cate, num = num)
          if(length(x1) > 1){
            x1 <- x1[[2]]
            supp1 <- .significance(x1)
          } else {
            supp1 <- 0 # It is the empty rule
          }
          
          if( supp1 >= mejor_soporte ){
            
            c1 <- .confianza(x1)
            c2 <-  .confianza(x)
            
            if( (supp1 > mejor_soporte) &&  c1 >= c2 ){
              mejor_soporte <- supp1
              mejor_regla <- regla_m
              mejor = TRUE
             }
          }
        }
        
      }
      
    }
    
  x1 <- .fit13(regla = mejor_regla, dataset = dataset, noClass = matrix(unlist(.separar(dataset)), nrow = length(dataset[[2]]) - 1, ncol = length(dataset[[7]])), targetClass = att_obj, por_cubrir = por_cubrir, n_Vars = dataset$nVars, nLabels = nLabels, max_regla = max_regla, marcar = TRUE, Objetivos = Objetivos, DNFRules = DNF_Rules, difuso = Objetivos[[4]], test = TRUE, cate = cate, num = num)[[2]]
  
  if(.confianza(x1) >= confianza_minima){
    mejor_regla
  } else {
    regla
  }
  
}

