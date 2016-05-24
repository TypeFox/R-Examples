
# This variable is for fix precision errors in numeric comparisons using the operator '==', '<=' or '>='
.tolerance <- 1e-10

#
#
# This function creates all the fuzzy intervals for a dataset
#
#
.create_fuzzyIntervals <- function(min, max, num_sets, types){
  #mapply(.FuzzyIntervals, min, max, num_sets, types, SIMPLIFY = FALSE)
  n_mats <- length(min)
  lst <- lapply(X = 1:n_mats, FUN = function(x) .FuzzyIntervals(min = min[x], max = max[x], num_sets = num_sets, types = types[x]))
  
  arr <- do.call(cbind, lst)
  dim(arr) <- c(num_sets, 3, n_mats)
  arr
}


#
#
#  Function to create fuzzy intervals over a single variable.
#
#

.FuzzyIntervals <- function(min, max, num_sets, types){
  
  xmin <- min - ( (max - min) / (num_sets - 1) )
  xmedio <- min
  xmax <- 0
  
  
  
  cjto_fuzzy <- matrix(nrow = num_sets, ncol = 3)
  
  for(j in 1:num_sets){ # For each fuzzy set
    if(types != "c"){ 
      xmax <- min + ( (max - min) / (num_sets - 1) ) * j
      
     
      
      #Save fuzzy set
      cjto_fuzzy[j, 1] <- xmin
      cjto_fuzzy[j, 2] <- xmedio
      cjto_fuzzy[j, 3] <- xmax
      
      #Modify variables 
      xmin <- xmedio
      xmedio <- xmax
    }
  }
  #Return
  cjto_fuzzy
  
}





#-----------------------------------------------------------------------------------
#    Creates crisp sets relative to the fuzzy sets generated before.
#    This works because the fuzzy sets are all triangular. If it does not, it does not work
#-----------------------------------------------------------------------------------

.createCrispIntervals <- function(fuzzyIntervals){
  n_mat <- dim(fuzzyIntervals)[3]
  n_vars <- dim(fuzzyIntervals)[1]
  crispIntervals <- lapply(X = 1:n_mat, FUN = function(x, fuzzy) .CrispIntervals(fuzzy[,,x]) , fuzzyIntervals)
  
  arr <- do.call(cbind, crispIntervals)
  dim(arr) <- c(n_vars, 2, n_mat)
  arr
}






#
#
# Create crisp sets for a single variable 
# 
#
.CrispIntervals <- function(fuzzyInterval){
  n_vars <- dim(fuzzyInterval)[1] # Cada fila de la matriz es una variable
  
  crispMatrix <- matrix(nrow = n_vars, ncol = 2)
  min <- fuzzyInterval[1,1]
  max <- (fuzzyInterval[1,3] + fuzzyInterval[1,2]) / 2
  crispMatrix[1,1] <- min
  crispMatrix[1,2] <- max
  for(i in 2:n_vars){
    min <- max
    max <- (fuzzyInterval[i,3] + fuzzyInterval[i,2]) / 2
    crispMatrix[i,1] <- min
    crispMatrix[i,2] <- max
  }
  crispMatrix
}





#
#
# Compute the belonging degree of an entire example set for the specified
#   fuzzy sets, this sets must be specified by means of their min, max 
#   and half value.
#   If you want to get the belonging degree for a subset of variables. Then
#   you must get only the values of the example corresponding to those variables
#   before execute this function.
#
#
.grado_pertenencia5 <- function(x, xmin, xmedio, xmax, n_matrices){
  
  x <- as.numeric(x)
  
  resultado <- numeric(length(x)) + 1
  xminX <- numeric(length(x)) + xmin + .tolerance
  xmedioX <- numeric(length(x)) + xmedio + .tolerance
  xmaxX <- numeric(length(x)) + xmax + .tolerance
  
  
  fuera_limites <- which( x <= xminX | x >= xmaxX )
  menorXMedio <- which( x < xmedioX & x > xminX )
  mayorXMedio <- which( x > xmedioX & x < xmaxX )
  
  resultado[fuera_limites] <- 0
  resultado[menorXMedio] <- (( x[menorXMedio] - xminX[menorXMedio] ) * (1 / (xmedioX[menorXMedio] - xminX[menorXMedio] )))
  resultado[mayorXMedio] <- (( xmaxX[mayorXMedio] - x[mayorXMedio] ) * (1 / (xmaxX[mayorXMedio] - xmedioX[mayorXMedio] )))
  
  resultado <- matrix(resultado, ncol = n_matrices, byrow = TRUE)
  
  resultado
}




#
#
# It is like grado_pertenencia5 but for crisp belonging and you must specify the min and max
#   value for crisp sets
#
#

.gradoPertenenciaCrisp2 <- function(x, xmin, xmax, DNF = FALSE){
  
  x <- as.numeric(x)
  result <- numeric(length(x))
  
 
    resulta <- which((x > xmin + .tolerance) & x <= (xmax + .tolerance))
    result[resulta] <- 1
    result <- matrix(result, ncol = length(xmax), byrow = TRUE)
    result
 
    
}





#
#
#  Return the maximum fuzzy degree for all variables that participe in a DNF rule.
#
#
.getMaxFuzzyForAVariable2 <- function(values, ejemplo_num){

  resultado <- .grado_pertenencia5(x = as.vector(ejemplo_num), xmin = values[2,], xmedio = values[3,], xmax = values[4,], n_matrices = NCOL(values))
 
  #Obtengo correctamente los grados de pertenencia, pero no consigo obtener el maximo de cada variable.  
  a <- which(!duplicated(values[1,]))
  long <- length(a)
  rangos <- vector(mode = "list", length = long)
  
  for(i in seq_len(long)){
    if(! is.na(a[i +1 ])){
      rangos[[i]]  <- a[i]:(a[i + 1] - 1)  
    } else {
      rangos[[i]] <- a[i]:NCOL(resultado)
    }
  }
  

if(NCOL(resultado) > 1){
  resultado <- t( apply(X = resultado, MARGIN = 1, FUN = function(x, rangos)  vapply(X = rangos, FUN = function(x, vector) max(vector[x]), FUN.VALUE = 1, x) , rangos) )
  resultado <- apply(X = resultado, MARGIN = 1, FUN = min)
} else {
  resultado <- t( apply(X = resultado, MARGIN = 1, FUN = function(x, rangos)  vapply(X = rangos, FUN = function(x, vector) max(vector[x]), FUN.VALUE = 1, x) , rangos) )
  
}
  resultado
  
}








#
#
# The same but for crips set on DNF Rules
#
#
.getMaxCrispForAVariable2 <- function(values, ejemplo_num){
  resultado <- .gradoPertenenciaCrisp2(x = as.vector(ejemplo_num), xmin = values[2,], xmax = values[3,], DNF = TRUE)
  #Obtengo correctamente los grados de pertenencia, pero no consigo obtener el maximo de cada variable.  
  a <- which(!duplicated(values[1,]))
  long <- length(a)
  rangos <- vector(mode = "list", length = long)
  
  for(i in seq_len(long)){
    if(! is.na(a[i +1 ])){
      rangos[[i]]  <- a[i]:(a[i + 1] - 1)  
    } else {
      rangos[[i]] <- a[i]:NCOL(resultado)
    }
  }
  
  if(NCOL(resultado) > 1){
    resultado <- t( apply(X = resultado, MARGIN = 1, FUN = function(x, rangos)  vapply(X = rangos, FUN = function(x, vector) max(vector[x]), FUN.VALUE = 1, x) , rangos) )
    resultado <- apply(X = resultado, MARGIN = 1, FUN = min)
  } else {
    resultado <- t( apply(X = resultado, MARGIN = 1, FUN = function(x, rangos)  vapply(X = rangos, FUN = function(x, vector) max(vector[x]), FUN.VALUE = 1, x) , rangos) )
    
  }
  resultado
  
}