#
#
#   THIS FILE CONTAINS SOME FUNCTIONS FOR DOING SOME STUFF INSIDE THE PACKAGE.
#   THIS MEANS THAT THE FUNCTIONS USED HERE MUST NOT BE EXPORTED. (EXCEPT ONE)
#
#


# Utils relative to obtain the fuzzy belonging degree

#--------------------------------------------------------------------



#'
#' Return the compatibility degrees of a rule with all instances of a given dataset.
#' 
#' The rules passed to this functions MUST have a vector representation in CANONICA form. This function
#' was made for being used mainly for the FuGePSD algorithm. 
#' 
#' @param ejemplo The instances of the dataset, a matrix with one example PER COLUMN and without the CLASS ATTRIBUTE.
#' @param rule_cat Part of the rule with the categorical values.
#' @param rule_num Part of the rule with the numerical values.
#' @param catParticip vector indicating which categorical attributes participate in the rule and so, they must be evaluated.
#' @param numParticip vector indicating which numerical attributes participate in the rule and so, they must be evaluated.
#' @param xmin numeric vector which indicate the minimum value of the fuzzy sets of every numeric attribute that participate in the rule.
#' @param xmedio numeric vector which indicate the medium value of the fuzzy sets of every numeric attribute that participate in the rule.
#' @param xmax numeric vector which indicate the maximum value of the fuzzy sets of every numeric attribute that participate in the rule.
#' @param n_matrices number of fuzzy sets that there are in the rule (The length of vectors xmin, xmedio and xmax)
#' @param max_cat numeric vector indicating the maximum value of categorical values.
#' @param max_num numeric vector indicating the maximum value for fuzzy partitions on every attribute
#' @param t_norm The T-norm to use to compute the compatibility degree. 0 for minimum t-norm. Other value for product t-norm
#' 
#' @return A vector with length 'number of examples' indicating their compatibility degree.
#' 
Rule.compatibility <- function(ejemplo, rule_cat, rule_num, catParticip, numParticip, xmin, xmedio, xmax, n_matrices,max_cat, max_num, t_norm){
  dispFuzzy <- numeric(NCOL(ejemplo)) + 1
  
  #Computation of membership degree.
  
  #Categorical variables.
  if(length(catParticip > 0)){
    ej_cat <- as.integer( ejemplo[catParticip,] )
    values <-  ceiling( ( which(ej_cat != rule_cat & ! (ej_cat == max_cat + 1 )) / (length(catParticip)) ) ) 
    #Examples no compatibles
    dispFuzzy[values] <- 0L
  }
  
  
  #Numerical Values 
  if(length(numParticip) > 0){
    ej_num <- as.vector( ejemplo[numParticip, which(dispFuzzy > 0) ] )
    
    #Fuzzy computation
    #Computes compatibility degree of every value of the whole dataset with the rule
    pertenencia <- .grado_pertenencia5(x = ej_num, xmin = xmin, xmedio = xmedio, xmax = xmax, n_matrices = n_matrices)
    
    if(t_norm == 0) { # MINIMUM T-NORM 
      dispFuzzy[which(dispFuzzy > 0)] <- apply(X = pertenencia, MARGIN = 1, FUN = min)
    } else { # PRODUCT T-NORM
      dispFuzzy[which(dispFuzzy > 0)] <- apply(X = pertenencia, MARGIN = 1, FUN = prod)
    }
  }
  
  
  dispFuzzy
  
}





#
#
#   Gets the membership degree of all examples in the dataset over a single rule. (Use with lapply)
#
#     ONLY FOR CAN RULES
# ejemplo es una matriz despues de usar .separar
# rule_cat is the categorcal variables that participate in the rule
# rule_num is the numerical variables that participate in the rule
# catParticip and numParticip are logical vectors for tell the function which rules of each type participe in the rule
# xmin, xmax, xmedio and xminCrisp and xmaxCrisp are the vectors with the fuzzy and crisp definition of the variables that participa in the rule
# max_cat is a vector the maximum value for categorical values and max_num is the same but for numerical variables.

.compara_CAN9 <- function(ejemplo, rule_cat, rule_num, catParticip, numParticip, xmin, xmedio, xmax, n_matrices, xminCrisp, xmaxCrisp, max_cat){
  dispFuzzy <- numeric(NCOL(ejemplo)) + 1
  dispCrisp <- integer(NCOL(ejemplo)) + 1L
  
  
  
  
  #Computation of membership degree.
  
  #Categorical variables.
  if(length(catParticip > 0)){
    ej_cat <- as.integer( ejemplo[catParticip,] )
    values <-  ceiling( ( which(ej_cat != rule_cat & ! (ej_cat == max_cat + 1 )) / (length(catParticip)) ) ) 
    #Examples no compatibles
    dispFuzzy[values] <- 0L
    dispCrisp[values] <- 0L
  }
  
  
  #Numerical Values 
  if(length(numParticip) > 0){
    ej_num <- as.vector( ejemplo[numParticip, which(dispFuzzy > 0) ] )
    #Fuzzy computation
    pertenencia <- .grado_pertenencia5(x = ej_num, xmin = xmin, xmedio = xmedio, xmax = xmax, n_matrices = n_matrices)
    dispFuzzy[which(dispFuzzy > 0)] <- apply(X = pertenencia, MARGIN = 1, FUN = min)
    
    #Crisp Computation
    pertenencia <- .gradoPertenenciaCrisp2(x = ej_num, xmin = xminCrisp, xmax = xmaxCrisp)
    dispCrisp[which(dispCrisp > 0)] <- apply(X = pertenencia, MARGIN = 1, FUN = min)
    
  }
  
  
  return(list( fuzzy = dispFuzzy, crisp = dispCrisp) )
  
}

#
#
#
# It works similar to compara_CAN9 but for DNF rules.
#

.comparaDNF4 <- function(ejemplo,  regla, regla_num, cat_particip, num_particip, max_regla_cat, max_regla_num, nLabels, fuzzySets, crispSet, valuesFuzzy, valuesCrisp){
  ejemplo_cat <- as.vector( ejemplo[cat_particip, ] )
  
  dispFuzzy <- numeric(NCOL(ejemplo)) + 1
  dispCrisp <- numeric(NCOL(ejemplo)) + 1
  #No salen los mismo resultados a?n.
  if(length(ejemplo_cat > 0)){  
    valCat <- (max_regla_cat + 1) + ejemplo_cat
    
    #Categorical Values
    fuera <- unique( ceiling(which(regla[valCat] == 0) / length(max_regla_cat)) )
    
    dispFuzzy[fuera] <- 0
    dispCrisp[fuera] <- 0
  }
  
  ejemplo_num <- ejemplo[num_particip, which(dispFuzzy > 0), drop = F] 
  
  #Numerical Values
  if(length(ejemplo_num) > 0){
    
    ejemplo_num <- ejemplo_num[valuesFuzzy[1,], ]
    #Fuzzy Computation
    dispFuzzy[which(dispFuzzy > 0)] <- .getMaxFuzzyForAVariable2(values = valuesFuzzy, ejemplo_num = ejemplo_num)
    #Crisp Computation
    dispCrisp[which(dispCrisp > 0)] <- .getMaxCrispForAVariable2(valuesCrisp, ejemplo_num)
  }
  
  list(fuzzy = dispFuzzy, crisp = dispCrisp)
  
}


#---------------------------------------------------------------------------
#   RETURN THE VALUES FOR CALCULATE THE QUALITY MEASURES
#
# - Return:
# -  [[1]] n(cond)  -> Ejemplos cubiertos por la regla
# -  [[2]] n(Tv ? cond) -> ejemplos cubiertos que cumplen el consecuente
# -  [[3]] FP -> Ejemplos que cumplen el antecedente pero no el consecuente (Falsos positivos)
# -  [[4]] Ns -> Numero de ejemplos en el dataset
# -  [[5]] n(TargetValue) -> numero de ejemplos que cumplen el consecuente
# -  [[6]] numero de ejemplos cubiertos de cada clase
# -  [[7]] numero de ejemplos de cada clase
# -  [[8]] ejemplos correctamente cubiertos que son nuevos
# -  [[9]] ejemplos de la clase objetivo que quedan por cubrir
# -  [[10]] suma difusa de los ejemplos cubiertos
# -  [[11]] suma difusa de los ejemplos correctamente cubiertos
# -  [[12]] suma difusa de los ejemplos nuevos correctamente cubiertos
# 
# ---------------------------------------------------------------


.get_values6 <- function(gr_perts, nombre_clases, dataset, targetClass, examples_perClass, cov, Ns, N_vars , por_cubrir, marcar = FALSE, test = FALSE, difuso = FALSE, NMEEF = FALSE){
  #Esto no es lo mejor, habr?a que buscar otra manera de utilizar la lista directamente
  dataset <- matrix(unlist(dataset), nrow = length(dataset[[1]]), ncol = length(dataset)) #MUCHO TIEMPO
  
  ejemplo_Cubiertos <- 0L
  sumaFuzzyejCubiertos <- 0
  ejemplosCorr_cubiertos <- 0L
  sumaFuzzyejCorrectamenteCubiertos <- 0
  ejemplosNuevos_Cubiertos <- 0L
  sumaFuzzyEjNuevos <- 0
  
  #Datos sobre el dataset
  
  #   cov_examplesFuzzy <- replicate(length(nombre_clases), 0)  # Para el calculo de la significancia
  #   names(cov_examplesFuzzy) <- nombre_clases
  
  cov_examplesCrisp <- integer(length(nombre_clases))  # Para el calculo de la significancia
  names(cov_examplesCrisp) <- nombre_clases
  
  fuzzyPerts <- gr_perts[[1]]
  crispPerts <- gr_perts[[2]]
  
  #Averiguamos los ejemplos cubiertos
  coveredFuzzy <- which( fuzzyPerts > 0)
  coveredCrisp <- which( crispPerts > 0)
  
  
  #Ejemplos cubiertor por la regla de cada clase (Significancia)
  #   tabla <- table( t( dataset[N_vars,coveredFuzzy]) )
  #   cov_examplesFuzzy[ names( tabla )] <- tabla 
  #tabla <- table( t( nombre_clases[ dataset[N_vars,coveredCrisp] + 1] ) )
  tabla <- improvedTable(dataset[, coveredCrisp, drop = F], nombre_clases)
  cov_examplesCrisp[names( tabla )] <- tabla 
  
  
  #Ejemplos cubiertos por la regla
  ejemplo_Cubiertos <- length(coveredCrisp)
  sumaFuzzyejCubiertos <- sum(fuzzyPerts[coveredFuzzy])
  
  
  #Ejemplos correctamente cubiertos
  p <- nombre_clases[ dataset[N_vars,coveredFuzzy] + 1] == targetClass 
  p1 <- nombre_clases[ dataset[N_vars,coveredCrisp] + 1] == targetClass 
  
  ejemplosCorr_cubiertos <- sum(p1)
  sumaFuzzyejCorrectamenteCubiertos <- sum(fuzzyPerts[coveredFuzzy[p]])
  
  
  #Ejemplos correctamente cubiertos que no estaban cubiertos anteriormente.  
  i <- cov[coveredFuzzy] == FALSE
  iC <- cov[coveredCrisp] == FALSE
  obj_notCoveredFuzzy <- which(p & i)
  obj_notCoveredCrisp <- which(p1 & iC)
  ejemplosNuevos_Cubiertos <- length(obj_notCoveredCrisp) #NCOL( dataset[ , obj_notCovered])
  sumaFuzzyEjNuevos <- sum(fuzzyPerts[coveredFuzzy[obj_notCoveredFuzzy]])
  
  #Marcar ejemplos nuevos cubiertos de la clase objetivo (Si es necesario)
  if(marcar){
    if(! difuso){
      cov[coveredCrisp[obj_notCoveredCrisp]] <- TRUE #Si se usa SOPORTE CRISP
    } else {
      cov[coveredFuzzy[p & i]] <- TRUE # USAR CUANDO SE USA SOPORTE DIFUSO
    }
    
    l <- list(ejemplo_Cubiertos, ejemplosCorr_cubiertos, NA, Ns, NROW(p[p]), cov_examplesCrisp, examples_perClass, ejemplosNuevos_Cubiertos, por_cubrir, sumaFuzzyejCubiertos, sumaFuzzyejCorrectamenteCubiertos, sumaFuzzyEjNuevos ) 
    conf <- .confianza(l)
    if( ! test) return(list(cov, conf)) 
    return(list(cov, l) )
  } else {
    
    #por_cubrir <- sum(obj_notCovered) 
    
    #Return 
    if(!NMEEF){
      return( list(ejemplo_Cubiertos, ejemplosCorr_cubiertos, NA, Ns, examples_perClass[[targetClass]], cov_examplesCrisp, examples_perClass, ejemplosNuevos_Cubiertos, por_cubrir, sumaFuzzyejCubiertos, sumaFuzzyejCorrectamenteCubiertos, sumaFuzzyEjNuevos ) )
    }else{
      cover <- (fuzzyPerts > 0 | crispPerts > 0) & nombre_clases[dataset[N_vars, ] + 1]== targetClass
      return( list(ejemplo_Cubiertos, ejemplosCorr_cubiertos, NA, Ns, examples_perClass[[targetClass]], cov_examplesCrisp, examples_perClass, ejemplosNuevos_Cubiertos, por_cubrir, sumaFuzzyejCubiertos, sumaFuzzyejCorrectamenteCubiertos, sumaFuzzyEjNuevos, cover ) )
    }
  }
  
}



.getVariableYValor <- function(value, max_valores){
  variable <- which( (value / max_valores) <= 1)[1]  # Solo queremos el primer valor. 
  vInicioVariable <- 1
  if(variable > 1){
    vInicioVariable <- max_valores[variable - 1] + 1
  }
  valor <- value - vInicioVariable + 1
  
  return(c(variable, valor))
}






















#---------------------------------------------------------------------
# oTHER UTILS

# 
# C.A.R. Hoare QuickSort Implementation
# @param v The vector to be ordered
# @param left First index of the subvector
# @param right Last index of the subvector
# @param index Index vector
# @return A list with two fields, vector which is the ordered vector and indices which is the sorted indexes of the original vector 
# 
.qsort <- function(v, left, right, index) {
 
  
  i = left
  j = right
  x = v[(left+right)/2]
  while(i <= j){
    while (v[i]<x && i<right)
      i <- i + 1
    while (x<v[j] && j>left)
      j <- j - 1
    if (i<=j) {
      y = v[i];
      v[i] = v[j];
      v[j] = y;
      aux = index[i];
      index[i] = index[j];
      index[j] = aux;
      i <- i + 1
      j <- j - 1
    }
}
  if (left<j){
   a <- .qsort(v,left,j,index)
    v[left:j] <- a$vector
    index[left:j] <- a$indices
  }
  if (i<right){
   b <- .qsort(v,i,right,index);
   v[i:right] <- b$vector
   index[i:right] <- b$indices
  }
 

  
  list(vector = v[left:right], indices = index[left:right])
}




#' Modifiy the number of Fuzzy Labels of the dataset.
#' 
#' This function change the number of fuzzy labels defined in the current KEEL dataset.
#' 
#' @details The fuzzy definitions used in the \code{keel} class are triangular.
#'     Because you can only specify the number of fuzzy definitions, all those definitions
#'     has the same width. With this function you can re-calculate this triangular fuzzy sets.
#' 
#' @param dataset The dataset to modify their fuzzy labels definitions. Must be a \code{keel} class.
#' @param nLabels The new number of fuzzy labels. An integer greater than zero.
#' 
#' @return  This function returns the same dataset with their fuzzy definitions modified.
#' 
#' @examples 
#' \dontrun{
#'     modifyFuzzyCrispIntervals(habermanTra, 2)
#'     modifyFuzzyCrispIntervals(habermanTra, 15)
#'}
#'

modifyFuzzyCrispIntervals <- function(dataset, nLabels){
    if(nLabels < 1)
      stop("The number of fuzzy sets ('nLabels') must be greater than zero.")
  
    dataset[["fuzzySets"]] <- .create_fuzzyIntervals(min = dataset$min, max = dataset$max, num_sets = nLabels, types = dataset$atributeTypes)
    dataset[["crispSets"]] <- .createCrispIntervals(fuzzyIntervals = dataset[["fuzzySets"]])
  
    dataset
}








#'
#' Change the target Variable of a \code{'keel'} Dataset
#' 
#' Change the actual target variable for another one if it is categorical.
#' 
#' @param dataset The KEEL dataset class
#' @param variable The position (or the name) of the variable to set as target Variable.
#' @return The dataset with the variables changed
#' 
#' 
#' @examples 
#' \dontrun{
#' changeTargetVariable(carTra, 3)
#' changeTargetVariable(carTra, "Doors")
#' 
#' Throws an error because the variable selected is numerical:
#' changeTargetVariable(habermanTra, 1)
#' }
#' 
#' 
changeTargetVariable <- function(dataset, variable){
  if(class(dataset) != "keel") stop( paste("'",substitute(dataset),"' is not a keel class", sep = ""))
  #if(variable >= dataset$nVars + 1) stop("variable is the same of the actual variable or is out of range")
  
  if(is.character(variable)){
    variable <- which(tolower(dataset$atributeNames) == tolower(variable))
    if(length(variable) == 0)
      stop(paste(variable, "is not a variable of this dataset."))
  }
  
  if(dataset[[3]][variable] != "c") stop("No categorical variable selected.")
  if(variable <= dataset$nVars){
  #Swap variables.
  dataset$data <- lapply(X = dataset$data , FUN = function(x, variable){ 
                       aux <- x[variable]; 
                       x[variable] <- x[length(x)]; 
                       x[length(x)] <- aux; 
                       x }, 
                       variable)
  
  #Swap Attribute Names
  aux <- dataset[[2]][variable]
  dataset[[2]][variable] <- dataset[[2]][length(dataset[[2]])]
  dataset[[2]][length(dataset[[2]])] <- aux
  
  #swap conjuntos
  dataset[["conjuntos"]][variable] <- dataset[["max"]][dataset[["nVars"]] + 1]
  
  #Swap Min
  aux <- dataset[[4]][variable]
  dataset[[4]][variable] <- dataset[[4]][length(dataset[[4]])]
  dataset[[4]][length(dataset[[4]])] <- aux
  
  #Swap Max
  aux <- dataset[[5]][variable]
  dataset[[5]][variable] <- dataset[[5]][length(dataset[[5]])]
  dataset[[5]][length(dataset[[5]])] <- aux
  
  #Change class_names Values
  dataset[["class_names"]] <- dataset[["categoricalValues"]][[variable]]
  
  #Swap categorical Values
  aux <- dataset[["categoricalValues"]][[variable]]
  dataset[["categoricalValues"]][[variable]] <- dataset[["categoricalValues"]][[length(dataset[["categoricalValues"]])]]
  dataset[["categoricalValues"]][[length(dataset[["categoricalValues"]])]] <- aux
  
  #Calculate new value for examplesPerClass
  clValues <- unlist(lapply(dataset$data, '[', dataset$nVars + 1))
  examplesPerClass <- lapply(X = seq_len(length(dataset$class_names)) - 1, FUN = function(x, data) sum(data == x), clValues)
  names(examplesPerClass) <- dataset$class_names
  dataset$examplesPerClass <- examplesPerClass
  
  
  }
  
  dataset
}







#
#
# Gets the variables that participate in a rule
#
#
.getParticipantes <- function(regla, max_regla, DNFRules){
  if(!DNFRules){
    participantes <- as.logical( (regla < max_regla) ) #Cuidado que ?sto puede que que no est? bien. No valor de la regla no es comparable con el ejemplo
  }else{
    
    participantes <- logical(length(max_regla) - 1)
    for(i in 2:length(max_regla)){
      ruleValues <- regla[(max_regla[i - 1] + 1):max_regla[i]]
      participantes[i-1] <- !(all(ruleValues == 1) | all(ruleValues == 0))
    }
  }
  
  participantes
  
}





#
# Returns de dataset without the last (class) column
#
.separar <- function(dataset){
  
  lapply(dataset$data, FUN = function(x) x[-length(x)])
  
}


#'
#' Returns the class attribute of the examples of a dataset.
#'
.getClassAttributes <- function(dataset){
  lapply(dataset, FUN = function(x) x[length(x)])
}



#
# returns de original dataes
#
.unir <- function(dataNoClass , classes){
  lapply(X = 1:length(dataNoClass), FUN = function(num, x,y) append(x[[num]],y[[num]]), dataNoClass, classes)
  
}









#
# Returns a matrix for select the variables that participate in a DNF rule for calculating their belonging degree
#
.getMatrixSelector <- function(regla_num, valor){
  
  values <- unlist(lapply(X = seq_len(length(regla_num)), 
                          FUN = function(x, regla, valor){
                            
                            a <- which(regla[[x]] > 0)
                            valores <- unlist(lapply(X = a, 
                                                     FUN = function(y, valor, mat){
                                                       c(y, valor, mat)
                                                     }, valor, x) )
                            valores
                          }, regla_num, valor))
  
  
  matrix(data = values, nrow = length(values) / 3, ncol = 3, byrow = TRUE)
  
}







#
#
# Obtain the fuzzy values of a DNF Rule
#
#
.getFuzzyValues <- function(regla_num, fuzzy,  crisp = FALSE){
  a <- .getMatrixSelector(regla_num = regla_num, valor = 1)
  variables <- a[,3]
  
  if(! crisp){
    xmin <- fuzzy[a]
    a[,2] <- 2
    xmedio <- fuzzy[a]
    a[,2] <- 3
    xmax <- fuzzy[a]
    
    rbind(variables, xmin, xmedio, xmax)
    
  } else {
    xmin <- fuzzy[a]
    a[,2] <- 2
    xmax <- fuzzy[a]
    
    rbind(variables, xmin, xmax)
  }
  
}






#
# Normalize a DNF by means of put all the non-participating variables filled with 0's.
#
.normalizeDNFRule <- function(regla, max){
  if(!anyNA(regla)){
    for(i in seq_len(length(max) - 1)){ 
      if(all(regla[(max[i] + 1):max[i+1]] == 1)){
        regla[(max[i] + 1):max[i+1]] <- 0
      }
    }
    
    regla
  } else{
    regla
  }
}















#' Launch a web interface for use the algorithms easily.
#' @description Launches a Shiny-based interface for the package in your browser.
#'     
#' @details The package \code{SDR} provide simple, shiny-based web interface for performs the taks 
#'     easily. The interface only work with new datasets loaded directly in the platform.
#'   
#'     The web application is structured as follows:
#' \itemize{
#'     \item{ The first you have to do is load your training and test files. This files must be valids KEEL format files.}
#'     \item{ After chose your datasets, you can view information about the dataset or execute the algorithm}
#'     \item{ You can choose the target variable or the variable to visualize and choose the target value or execute the algorithm for all the values.}
#'     \item{ Choosed the target variable, you can  choose the algorithm to execute and change his parameters with the controls provided.}
#'     \item{ After you can execute the algorithm. The results are exposed in three tabs that are at the top of the page, just at the right of the "Exploratory Analysis" tab.}
#' }
#'     The tables can be sorted for each value and also you can search and filter values.
#'     
#'     
#' @examples
#'\dontrun{
#' library(SDR)
#' SDR_GUI()
#'}
#'     
#' @export
 SDR_GUI <- function(){
   packages <- installed.packages()[,1]
   if(! "shiny" %in% packages){
     if(tolower(.yesno("Package 'shiny' is not installed and must be installed to run this GUI. Do you want to install it? (Y/n): ")) == "y"){
       install.packages("shiny")
       cat("Launching interface...")
       shiny::runApp(appDir = system.file("shiny", package="SDR"), launch.browser = TRUE)
       
       invisible()
     } else {
       cat("Package not installed. Execution aborted.")
     }
   } else {
     shiny::runApp(appDir = system.file("shiny", package="SDR"), launch.browser = TRUE)
     
     invisible()
   }
 }





# Make a yes/no question to the user
.yesno <- function(texto){
  line <- readline(texto)
  line
}















######################################################################
#                                                                    #
#                   MATHEMATICAL UTILS                               #    
#                                                                    #
#                                                                    #
######################################################################


#'
#' returns an integer between low and high. EXCLUDING high
#' @param low Lower bound (included)
#' @param high Upper bound (NOT included)
#' @return a uniform-distributed integer value in [low, high)
#' 
.randInt <- function(low, high){
  floor( low + (high - low) * runif(1) )
}



#'
#' returns an integer between low and high. INCLUDING high
#' @param low Lower bound (included)
#' @param high Upper bound (included)
#' @return a uniform-distributed integer value in [low, high)
#' 
#' 
.randIntClosed <- function(low, high){
  floor( low + ((high + 1) - low) * runif(1) )
}



#'
#' returns an integer between low and high. EXCLUDING low and high
#' @param low Lower bound (NOT included)
#' @param high Upper bound (NOT included)
#' @return a uniform-distributed integer value in [low, high)
#' 
.randIntClosed <- function(low, high){
  floor( (low+1) + (high - (low+1)) * runif(1) )
}


#'
#' returns a number between low and high. EXCLUDING high
#' @param low Lower bound (included)
#' @param high Upper bound (NOT included)
#' @return a uniform-distributed integer value in [low, high)
#' 
.randDouble <- function(low, high){
  low + (high - low) * runif(1) 
}



#'
#' returns a number between low and high. Including high and EXCLUDING excluded
#' @param low Lower bound (included)
#' @param high Upper bound (NOT included)
#' @param excluded. The number to exclude, it does not check if it is in the range [low,high]
#' @return a uniform-distributed integer value in [low, high)
#' 
.randIntExcluded <- function(low, high, excluded){
  number <- .randIntClosed(low, high)
  while(number == excluded){
    number <- .randIntClosed(low, high)
  }
  number
}




#'
#' Parse a time differente to "x hours, y minutes and z seconds"
#'
#' @param actual End time in UNIX int format (i.e. as.numeric(Sys.time()))
#' @param initial Initial time in Unix format.
#'
#' @return A human-readable string with time difference.
#'
parseTime <- function(actual, initial){
  dif <- actual - initial
  horas <- 0
  minutos <- 0
  segundos <- 0
  
  if(dif >= 3600){
    horas <- floor(dif / 3600)
    dif <- dif %% 3600
  }
  
  if(dif >= 60){
    minutos <- floor(dif / 60)
    segundos <- dif %% 60
  } else {
    segundos <- dif
  }
  
  
  paste(horas, " hours, ", minutos, " minutes and ", round(segundos, 2) , " seconds.", sep = "")
}


#'
#' Improved table creation for .get_values6
#' 
#' @param dataset A matrix with the data
#' @param classNames a vector with the names of the attributes.
#' 
#' @return a named vector with the number of instances per class.
#' 
improvedTable <- function(dataset, classNames){
  tabla <-
    vapply(
      X = seq_len(length(classNames)) - 1, FUN = function(x, data)
        sum(data == x), integer(1), dataset[nrow(dataset),]
    )
  names(tabla) <- classNames
  tabla
}



#' S3 function to summary a keel object
#' 
#' Summary relevant data of a \code{keel} dataset.
#' 
#' @param object A \code{keel} class.
#' @param ... Additional arguments to the summary function.
#' 
#' @details This function show important information about the \code{keel} dataset for the user. Note that it does not 
#' show all the information available. The rest is only for the algorithms. The values that appear are accesible by the
#' \code{$} operator, e.g. dataset$relation or dataset$examplesPerClass.
#' 
#'@examples 
#'  
#'summary(carTra) 
#' 
#'@export
summary.keel <- function(object, ...){
  cat(paste("Summary of the keel object: '", substitute(object),"'", sep = ""),
      paste("\t- relation:", object$relation),
      paste("\t- nVars:", object$nVars),
      paste("\t- Ns:", object$Ns),
      paste("\t- attributeNames:", paste(object$atributeNames, collapse = ", ")),
      paste("\t- class_names:", paste(object$class_names, collapse = ", ")),
      paste("\t- examplesPerClass:" ,paste(unlist(object$examplesPerClass), collapse = ", "))
      , sep = "\n")
}



#'  S3 function to print in console the contents of the dataset
#'  
#'  This function shows the matrix of data uncoded.
#'  
#'  @param x The \code{keel} object to view
#'  @param ... Additional arguments passed to the print function
#'  
#'  @details This function show the matix of data. Internally, a \code{keel} object has a list of of examples
#'  and this examples are coded numerically. This function decode these examples and convert the list into a matrix.
#'  
#'  @return a matrix with the dataset uncoded.
#'  
#' @export
print.keel <- function(x, ...){
  data <- lapply(x$data,
                 function(x, categoricos)
                   vapply(seq_len(length(x)), function(i, ejemplo, cateValues){
                     if(is.na(cateValues[[i]][1])){
                       as.character(ejemplo[i])
                     } else{
                       cateValues[[i]][ejemplo[i] + 1]
                     }
                   }, character(1), x, categoricos)
                 
                 , x$categoricalValues
                 
  )
  
  print(matrix(data = unlist(data), ncol = x$nVars + 1, byrow = TRUE, dimnames = list(NULL,x$atributeNames)), ...)
}
