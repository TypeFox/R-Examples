#' @name byte_mutation
#' @title Performs mutation operation on a given double vector
#' @description This function is not called directly but is given as a parameter in \code{GA::ga} function. 
#' In \code{GA::ga}, if the parameter \code{mutation=} is set to \code{byte_mutation} than
#' the byte-coded mutation operator is applied in the genetic search. In \code{mcga2} function, the hard-coded 
#' mutation parameter is set to byte_mutation by definition. Byte-mutation function simply takes an double vector and
#' changes bytes of this values by +1 or -1 using the pre-determined mutation probabilty. 
#' @param object A \code{GA::ga} object
#' @param parent Index of the candidate solution of the current population
#' @param ... Additional arguments to be passed to the function
#' @return Mutated double vector
#' @author Mehmet Hakan Satman - mhsatman@istanbul.edu.tr
#' @references M.H.Satman (2013), Machine Coded Genetic Algorithms for Real Parameter Optimization Problems, Gazi University Journal of Science, Vol 26, No 1, pp. 85-95
#' @examples
#' f <- function(x){ 
#'   return(-sum( (x-5)^2 ) )
#' }
#' myga <- GA::ga(type="real-valued", fitness = f, popSize = 100, maxiter = 200, 
#'               min = rep(-50,5), max = rep(50,5), crossover = byte_crossover,
#'               mutation = byte_mutation)
#' print(myga@solution)
byte_mutation <- function (object, parent, ...){
  mutate <- as.vector(object@population[parent, ])
  mutate <- ByteCodeMutationUsingDoubles(mutate, object@pmutation)
  EnsureBounds(mutate,object@min, object@max)
  return(mutate)
}

#' @name byte_mutation_dynamic
#' @title Performs mutation operation on a given double vector using dynamic mutation probabilities
#' @description This function is not called directly but is given as a parameter in \code{GA::ga} function. 
#' In \code{GA::ga}, if the parameter \code{mutation=} is set to \code{byte_mutation_dynamic} than
#' the byte-coded mutation operator is applied in the genetic search. In \code{mcga2} function, the hard-coded 
#' mutation parameter is set to byte_mutation by definition. Byte-mutation function simply takes an double vector and
#' changes bytes of this values by +1 or -1 using the dynamically decreased and pre-determined mutation probabilty. 
#' @param object A \code{GA::ga} object
#' @param parent Index of the candidate solution of the current population
#' @param ... Additional arguments to be passed to the function
#' @return Mutated double vector
#' @author Mehmet Hakan Satman - mhsatman@istanbul.edu.tr
#' @references M.H.Satman (2013), Machine Coded Genetic Algorithms for Real Parameter Optimization Problems, Gazi University Journal of Science, Vol 26, No 1, pp. 85-95
#' @examples
#' f <- function(x){ 
#'   return(-sum( (x-5)^2 ) )
#' }
#' myga <- GA::ga(type="real-valued", fitness = f, popSize = 100, maxiter = 200, 
#'               min = rep(-50,5), max = rep(50,5), crossover = byte_crossover,
#'               mutation = byte_mutation_dynamic, pmutation = 0.10)
#' print(myga@solution)
byte_mutation_dynamic <- function (object, parent, ...){
  mutate <- as.vector(object@population[parent, ])
  pmutation <- object@pmutation - object@pmutation * (object@iter / object@maxiter) 
  mutate <- ByteCodeMutationUsingDoubles(mutate, pmutation)
  EnsureBounds(mutate,object@min, object@max)
  return(mutate)
}


#' @name byte_mutation_random
#' @title Performs mutation operation on a given double vector
#' @description This function is not called directly but is given as a parameter in \code{GA::ga} function. 
#' In \code{GA::ga}, if the parameter \code{mutation=} is set to \code{byte_mutation_random} than
#' the byte-coded mutation operator is applied in the genetic search. In \code{mcga2} function, the hard-coded 
#' mutation parameter is set to byte_mutation by definition. This function simply takes an double vector and
#' changes bytes randomly in the range of [0,255] using the pre-determined mutation probabilty. 
#' @param object A \code{GA::ga} object
#' @param parent Index of the candidate solution of the current population
#' @param ... Additional arguments to be passed to the function
#' @return Mutated double vector
#' @author Mehmet Hakan Satman - mhsatman@istanbul.edu.tr
#' @references M.H.Satman (2013), Machine Coded Genetic Algorithms for Real Parameter Optimization Problems, Gazi University Journal of Science, Vol 26, No 1, pp. 85-95
#' @examples
#' f <- function(x){ 
#'   return(-sum( (x-5)^2 ) )
#' }
#' myga <- GA::ga(type="real-valued", fitness = f, popSize = 100, maxiter = 200, 
#'               min = rep(-50,5), max = rep(50,5), crossover = byte_crossover,
#'               mutation = byte_mutation_random, pmutation = 0.20)
#' print(myga@solution)
byte_mutation_random <- function (object, parent, ...){
  mutate <- as.vector(object@population[parent, ])
  mutate <- ByteCodeMutationUsingDoublesRandom(mutate, object@pmutation)
  EnsureBounds(mutate,object@min, object@max)
  return(mutate)
}



#' @name byte_mutation_random_dynamic
#' @title Performs mutation operation on a given double vector with dynamic mutation probabilities
#' @description This function is not called directly but is given as a parameter in \code{GA::ga} function. 
#' In \code{GA::ga}, if the parameter \code{mutation=} is set to \code{byte_mutation_random_dynamic} than
#' the byte-coded mutation operator with dynamic probabilities is applied in the genetic search. In \code{mcga2} function, the hard-coded 
#' mutation parameter is set to byte_mutation by definition. This function simply takes an double vector and
#' changes bytes randomly in the range of [0,255] using the descrasing values of pre-determined mutation probabilty by generations. 
#' @param object A \code{GA::ga} object
#' @param parent Index of the candidate solution of the current population
#' @param ... Additional arguments to be passed to the function
#' @return Mutated double vector
#' @author Mehmet Hakan Satman - mhsatman@istanbul.edu.tr
#' @references M.H.Satman (2013), Machine Coded Genetic Algorithms for Real Parameter Optimization Problems, Gazi University Journal of Science, Vol 26, No 1, pp. 85-95
#' @examples
#' f <- function(x){ 
#'   return(-sum( (x-5)^2 ) )
#' }
#' # Increase popSize and maxiter for more precise solutions
#' myga <- GA::ga(type="real-valued", fitness = f, popSize = 100, maxiter = 200, 
#'               min = rep(-50,5), max = rep(50,5), crossover = byte_crossover,
#'               mutation = byte_mutation_random_dynamic, pmutation = 0.20)
#' print(myga@solution)
byte_mutation_random_dynamic <- function (object, parent, ...){
  mutate <- as.vector(object@population[parent, ])
  pmutation <- object@pmutation - object@pmutation * (object@iter / object@maxiter) 
  mutate <- ByteCodeMutationUsingDoublesRandom(mutate, pmutation)
  EnsureBounds(mutate,object@min, object@max)
  return(mutate)
}

#' @name byte_crossover
#' @title Performs crossover operation on a pair of two selected parent candidate solutions
#' @description This function is not called directly but is given as a parameter in \code{GA::ga} function. 
#' In \code{GA::ga}, if the parameter \code{crossover=} is set to \code{byte_crossover} than
#' the byte-coded crossover operator is applied in the genetic search. In \code{mcga2} function, the hard-coded 
#' crossover parameter is set to byte_crossover by definition. \code{byte_crossover} function simply takes two double vectors 
#' (parents) and combines the bytes of doubles using a Uniform distribution with parameters 0 and 1. 
#' @param object A \code{GA::ga} object
#' @param parents Indices of the selected parents 
#' @param ... Additional arguments to be passed to the function
#' @return List of two generated offspring
#' @author Mehmet Hakan Satman - mhsatman@istanbul.edu.tr
#' @references M.H.Satman (2013), Machine Coded Genetic Algorithms for Real Parameter Optimization Problems, Gazi University Journal of Science, Vol 26, No 1, pp. 85-95
#' @seealso mcga2
#' @examples
#' f <- function(x){ 
#'   return(-sum( (x-5)^2 ) )
#' }
#' myga <- GA::ga(type="real-valued", fitness = f, popSize = 100, maxiter = 200, 
#'               min = rep(-50,5), max = rep(50,5), crossover = byte_crossover,
#'               mutation = byte_mutation)
#' print(myga@solution)
byte_crossover <- function (object, parents, ...){
  parents <- object@population[parents, , drop = FALSE]
  ch1 <- parents[1,]
  ch2 <- parents[2,]
  
  offs <- UniformCrossOverOnDoublesUsingBytes(ch1, ch2)
  off1 <- offs[[1]]
  off2 <- offs[[2]]
  EnsureBounds(off1, object@min, object@max)
  EnsureBounds(off2, object@min, object@max)
  
  out <- list(children = rbind(off1,off2), fitness = rep(NA, 2))
  return(out)
}



#' @name byte_crossover_1p
#' @title Performs one-point crossover operation on a pair of two selected parent candidate solutions
#' @description This function is not called directly but is given as a parameter in \code{GA::ga} function. 
#' In \code{GA::ga}, if the parameter \code{crossover=} is set to \code{byte_crossover_1p} than
#' the byte-coded one-point crossover operator is applied in the genetic search. In \code{mcga2} function, the hard-coded 
#' crossover parameter is set to byte_crossover by definition. \code{byte_crossover_1p} function simply takes two double vectors 
#' (parents) and combines the bytes of doubles using given cut-point. 
#' @param object A \code{GA::ga} object
#' @param parents Indices of the selected parents 
#' @param ... Additional arguments to be passed to the function
#' @return List of two generated offspring
#' @author Mehmet Hakan Satman - mhsatman@istanbul.edu.tr
#' @references M.H.Satman (2013), Machine Coded Genetic Algorithms for Real Parameter Optimization Problems, Gazi University Journal of Science, Vol 26, No 1, pp. 85-95
#' @seealso mcga2
#' @examples
#' f <- function(x){ 
#'   return(-sum( (x-5)^2 ) )
#' }
#' myga <- GA::ga(type="real-valued", fitness = f, popSize = 100, maxiter = 200, 
#'               min = rep(-50,5), max = rep(50,5), crossover = byte_crossover_1p,
#'               mutation = byte_mutation)
#' print(myga@solution)
byte_crossover_1p <- function (object, parents, ...){
  parents <- object@population[parents, , drop = FALSE]
  ch1 <- parents[1,]
  ch2 <- parents[2,]
  offs <- OnePointCrossOverOnDoublesUsingBytes(ch1, ch2, sample(1:(length(ch1) * SizeOfDouble()),1)[1])
  off1 <- offs[[1]]
  off2 <- offs[[2]]
  EnsureBounds(off1, object@min, object@max)
  EnsureBounds(off2, object@min, object@max)
  out <- list(children = rbind(off1,off2), fitness = rep(NA, 2))
  return(out)
}


#' @name byte_crossover_2p
#' @title Performs two-point crossover operation on a pair of two selected parent candidate solutions
#' @description This function is not called directly but is given as a parameter in \code{GA::ga} function. 
#' In \code{GA::ga}, if the parameter \code{crossover=} is set to \code{byte_crossover_2p} than
#' the byte-coded two-point crossover operator is applied in the genetic search. In \code{mcga2} function, the hard-coded 
#' crossover parameter is set to byte_crossover by definition. \code{byte_crossover_2p} function simply takes two double vectors 
#' (parents) and combines the bytes of doubles using given cutpoint1 and cutpoint2. 
#' @param object A \code{GA::ga} object
#' @param parents Indices of the selected parents 
#' @param ... Additional arguments to be passed to the function
#' @return List of two generated offspring
#' @author Mehmet Hakan Satman - mhsatman@istanbul.edu.tr
#' @references M.H.Satman (2013), Machine Coded Genetic Algorithms for Real Parameter Optimization Problems, Gazi University Journal of Science, Vol 26, No 1, pp. 85-95
#' @seealso mcga2
#' @examples
#' f <- function(x){ 
#'   return(-sum( (x-5)^2 ) )
#' }
#' myga <- GA::ga(type="real-valued", fitness = f, popSize = 100, maxiter = 200, 
#'               min = rep(-50,5), max = rep(50,5), crossover = byte_crossover_2p,
#'               mutation = byte_mutation)
#' print(myga@solution)
byte_crossover_2p <- function (object, parents, ...){
  parents <- object@population[parents, , drop = FALSE]
  ch1 <- parents[1,]
  ch2 <- parents[2,]
  cutpoints <- sort(sample(1:(length(ch1)*SizeOfDouble()), 2, replace = FALSE))
  offs <- TwoPointCrossOverOnDoublesUsingBytes(ch1, ch2, cutpoints[1], cutpoints[2])
  off1 <- offs[[1]]
  off2 <- offs[[2]]
  EnsureBounds(off1, object@min, object@max)
  EnsureBounds(off2, object@min, object@max)
  out <- list(children = rbind(off1,off2), fitness = rep(NA, 2))
  return(out)
}


#' @name sbx_crossover
#' @title Performs sbx (simulated binary) crossover operation on a pair of two selected parent candidate solutions
#' @description This function is not called directly but is given as a parameter in \code{GA::ga} function. 
#' In \code{GA::ga}, if the parameter \code{crossover=} is set to \code{sbx_crossover} than
#' the sbx crossover operator is applied in the genetic search. sbx_crossover mimics the classical single-point crossover operator
#' in binary genetic algorithms. 
#' @param object A \code{GA::ga} object
#' @param parents Indices of the selected parents 
#' @param ... Additional arguments to be passed to the function
#' @return List of two generated offspring
#' @author Mehmet Hakan Satman - mhsatman@istanbul.edu.tr
#' @references Deb, Kalyanmoy, and Ram Bhushan Agrawal. "Simulated binary crossover for continuous search space." Complex systems 9.2 (1995): 115-148.
#' @examples
#' f <- function(x){ 
#'   return(-sum( (x-5)^2 ) )
#' }
#' myga <- ga(type="real-valued", fitness = f, popSize = 100, maxiter = 100, 
#'            min = rep(-50,5), max = rep(50,5), crossover = sbx_crossover)
#' print(myga@solution)
sbx_crossover <- function (object, parents, ...) 
{
  parents <- object@population[parents, , drop = FALSE]
  p <- dim(parents)[2]
  if(!exists("nc")){
    nc <- 50
  }
  u <- runif(p)
  betaq <- rep(0,p)
  for (i in 1:p){
    if(u[i] <= 0.5){
      betaq[i] <- (2*u[i]) ^ (1/(nc + 1))
    }else{
      betaq[i] <- (1/(2*(1-u[i]))) ^ (1/(nc+1))
    }
  }
  off1 <- 0.5 * ((1+betaq)*parents[1,] + (1-betaq[i])*parents[2,])
  off2 <- 0.5 * ((1-betaq)*parents[1,] + (1+betaq[i])*parents[2,])
  out <- list(children = rbind(off1,off2), fitness = rep(NA, 2))
  return(out)
}

#' @name flat_crossover
#' @title Performs flat crossover operation on a pair of two selected parent candidate solutions
#' @description This function is not called directly but is given as a parameter in \code{GA::ga} function. 
#' In \code{GA::ga}, if the parameter \code{crossover=} is set to \code{flat_crossover} than
#' the flat crossover operator is applied in the genetic search. \code{flat_crossover} draws a random number between parents' genes and returns a pair of generated offspring
#' @param object A \code{GA::ga} object
#' @param parents Indices of the selected parents 
#' @param ... Additional arguments to be passed to the function
#' @return List of two generated offspring
#' @author Mehmet Hakan Satman - mhsatman@istanbul.edu.tr
#' @examples
#' f <- function(x){ 
#'   return(-sum( (x-5)^2 ) )
#' }
#' myga <- ga(type="real-valued", fitness = f, popSize = 100, maxiter = 100, 
#'            min = rep(-50,5), max = rep(50,5), crossover = flat_crossover)
#' print(myga@solution)
flat_crossover <- function (object, parents, ...){
  parents <- object@population[parents, , drop = FALSE]
  p <- dim(parents)[2]
  off1 <- rep(0,p)
  off2 <- rep(0,p)
  for (i in 1:p){
    my.min <- min(parents[1,i], parents[2,i])
    my.max <- max(parents[1,i], parents[2,i])
    off1[i] <- runif(1, min=my.min, max = my.max)
    off2[i] <- runif(1, min=my.min, max = my.max)
  }
  out <- list(children = rbind(off1,off2), fitness = rep(NA, 2))
  return(out)
}


#' @name arithmetic_crossover
#' @title Performs arithmetic crossover operation on a pair of two selected parent candidate solutions
#' @description This function is not called directly but is given as a parameter in \code{GA::ga} function. 
#' In \code{GA::ga}, if the parameter \code{crossover=} is set to \code{arithmetic_crossover} than
#' the arithmetic crossover operator is applied in the genetic search. \code{arithmetic_crossover} generates offspring using the weighted mean of parents' genes. Weights are drawn randomly.
#' @param object A \code{GA::ga} object
#' @param parents Indices of the selected parents 
#' @param ... Additional arguments to be passed to the function
#' @return List of two generated offspring
#' @author Mehmet Hakan Satman - mhsatman@istanbul.edu.tr
#' @examples
#' f <- function(x){ 
#'   return(-sum( (x-5)^2 ) )
#' }
#' myga <- ga(type="real-valued", fitness = f, popSize = 100, maxiter = 100, 
#'            min = rep(-50,5), max = rep(50,5), crossover = arithmetic_crossover)
#' print(myga@solution)
arithmetic_crossover <- function (object, parents, ...){
  parents <- object@population[parents, , drop = FALSE]
  p <- dim(parents)[2]
  alpha <- runif(1,min=0.0,max=1.0)
  off1 <- alpha * parents[1,] + (1-alpha) * parents[2,]
  off2 <- (1-alpha) * parents[1,] + (alpha) * parents[2,]
  out <- list(children = rbind(off1,off2), fitness = rep(NA, 2))
  return(out)
}

#' @name blx_crossover
#' @title Performs blx (blend) crossover operation on a pair of two selected parent candidate solutions
#' @description This function is not called directly but is given as a parameter in \code{GA::ga} function. 
#' In \code{GA::ga}, if the parameter \code{crossover=} is set to \code{blx_crossover} than
#' the blx crossover operator is applied in the genetic search. 
#' @param object A \code{GA::ga} object
#' @param parents Indices of the selected parents 
#' @param ... Additional arguments to be passed to the function
#' @return List of two generated offspring
#' @author Mehmet Hakan Satman - mhsatman@istanbul.edu.tr
#' @examples
#' f <- function(x){ 
#'   return(-sum( (x-5)^2 ) )
#' }
#' myga <- ga(type="real-valued", fitness = f, popSize = 100, maxiter = 100, 
#'            min = rep(-50,5), max = rep(50,5), crossover = blx_crossover)
#' print(myga@solution)
blx_crossover <- function (object, parents, ...){
  parents <- object@population[parents, , drop = FALSE]
  p <- dim(parents)[2]
  alpha <- 0.5
  off1 <- rep(0,p)
  off2 <- rep(0,p)
  for (i in 1:p){
    r.min <- min(parents[1,i], parents[2,i])
    r.max <- max(parents[1,i], parents[2,i])
    I <- r.max - r.min
    off1[i] <- runif(1, min=r.min - I*alpha, max=r.max + I*alpha)
    off2[i] <- runif(1, min=r.min - I*alpha, max=r.max + I*alpha)
  }
  out <- list(children = rbind(off1,off2), fitness = rep(NA, 2))
  return(out)
}

#' @name linear_crossover
#' @title Performs linear crossover operation on a pair of two selected parent candidate solutions
#' @description This function is not called directly but is given as a parameter in \code{GA::ga} function. 
#' In \code{GA::ga}, if the parameter \code{crossover=} is set to \code{linear_crossover} than
#' the linear crossover operator is applied in the genetic search. \code{linear_crossover} generates three offspring and performs a selection mechanism to determine best two of them. 
#' @param object A \code{GA::ga} object
#' @param parents Indices of the selected parents 
#' @param ... Additional arguments to be passed to the function
#' @return List of two generated offspring
#' @author Mehmet Hakan Satman - mhsatman@istanbul.edu.tr
#' @examples
#' f <- function(x){ 
#'   return(-sum( (x-5)^2 ) )
#' }
#' myga <- ga(type="real-valued", fitness = f, popSize = 100, maxiter = 100, 
#'            min = rep(-50,5), max = rep(50,5), crossover = linear_crossover)
#' print(myga@solution)
linear_crossover <- function (object, parents, ...){
  parents <- object@population[parents, , drop = FALSE]
  p <- dim(parents)[2]
  
  off1 <- (1/2) * parents[1,] + (1/2) * parents[2,]
  off2 <- (3/2) * parents[1,] - (1/2) * parents[2,]
  off3 <- (-1/2)* parents[1,] + (3/2) * parents[2,]
  
  all.off <-  list(off1,off2,off3)
  fit <- unlist(lapply(X = all.off, object@call$fitness))
  max.index <- which.max(fit)
  best1 <- all.off[[max.index]]
  best2 <- all.off[[setdiff(1:3, max.index)]]
  out <- list(children = rbind(best1,best2), fitness = rep(NA, 2))
  return(out)
}

#' @name unfair_average_crossover
#' @title Performs unfair average crossover operation on a pair of two selected parent candidate solutions
#' @description This function is not called directly but is given as a parameter in \code{GA::ga} function. 
#' In \code{GA::ga}, if the parameter \code{crossover=} is set to \code{unfair_average_crossover} than
#' the unfair average crossover operator is applied in the genetic search. 
#' @param object A \code{GA::ga} object
#' @param parents Indices of the selected parents 
#' @param ... Additional arguments to be passed to the function
#' @return List of two generated offspring
#' @author Mehmet Hakan Satman - mhsatman@istanbul.edu.tr
#' @examples
#' f <- function(x){ 
#'   return(-sum( (x-5)^2 ) )
#' }
#' myga <- ga(type="real-valued", fitness = f, popSize = 100, maxiter = 100, 
#'            min = rep(-50,5), max = rep(50,5), crossover = unfair_average_crossover)
#' print(myga@solution)
unfair_average_crossover <- function (object, parents, ...){
  parents <- object@population[parents, , drop = FALSE]
  p <- dim(parents)[2]
  off1 <- rep(0,p)
  off2 <- rep(0,p)
  alpha <- runif(1, min=0.0, max=0.5)
  j <- sample(1:p,1)
  for (i in 1:j){
    if (i <= j){
      off1[i] <- (1+alpha) * parents[1,i] - alpha * parents[2,i]
      off2[i] <- -alpha * parents[1,i] + (1+alpha) * parents[2,i]
    }else{
      off1[i] <- (1-alpha) * parents[1,i] + alpha * parents[2,i]
      off2[i] <- alpha * parents[1,i] + (1-alpha) * parents[2,i]
    }
  }
  out <- list(children = rbind(off1,off2), fitness = rep(NA, 2))
  return(out)
}



#' @name mcga2
#' @title Performs a machine-coded genetic algorithm search for a given optimization problem
#' @description \code{mcga2} is the improvement version of the standard mcga function as it is based on the \code{GA::ga} function. The 
#' \code{byte_crossover} and the \code{byte_mutation} operators are the main reproduction operators and these operators uses the byte 
#' representations of parents in the computer memory. 
#' @param fitness The goal function to be maximized
#' @param ... Additional arguments to be passed to the fitness function
#' @param min Vector of lower bounds of variables
#' @param max Vector of upper bounds of variables
#' @param population Initial population. It is \code{gaControl("real-valued")$population} by default.
#' @param selection Selection operator. It is \code{gaControl("real-valued")$selection} by default. 
#' @param crossover Crossover operator. It is \code{byte_crossover} by default.
#' @param mutation Mutation operator. It is \code{byte_mutation} by default. Other values can be given including \code{byte_mutation_random}, 
#' \code{byte_mutation_dynamic} and \code{byte_mutation_random_dynamic}
#' @param popSize Population size. It is 50 by default
#' @param pcrossover Probability of crossover. It is 0.8 by default
#' @param pmutation Probability of mutation. It is 0.1 by default
#' @param elitism Number of elitist solutions. It is \code{base::max(1, round(popSize*0.05))} by default
#' @param maxiter Maximum number of generations. It is 100 by default
#' @param run The genetic search is stopped if the best solution has not any improvements in last \code{run} generations. By default it is \code{maxiter}
#' @param maxfitness Upper bound of the fitness function. By default it is Inf
#' @param names Vector of names of the variables. By default it is \code{NULL}
#' @param parallel If TRUE, fitness calculations are performed parallel. It is FALSE by default
#' @param monitor The monitoring function for printing some information about the current state of the genetic search. It is \code{gaMonitor} by default
#' @param seed The seed for random number generating. It is \code{NULL} by default
#' @return Returns an object of class \code{ga-class}
#' @seealso GA::ga
#' @author Mehmet Hakan Satman - mhsatman@istanbul.edu.tr
#' @references M.H.Satman (2013), Machine Coded Genetic Algorithms for Real Parameter Optimization Problems, Gazi University Journal of Science, Vol 26, No 1, pp. 85-95
#' @references Luca Scrucca (2013). GA: A Package for Genetic Algorithms in R. Journal of Statistical Software, 53(4), 1-37. URL \url{http://www.jstatsoft.org/v53/i04/}
#' @examples
#' f <- function(x){ 
#'   return(-sum( (x-5)^2 ) )
#' }
#' myga <- mcga2(fitness = f, popSize = 100, maxiter = 300, 
#'               min = rep(-50,5), max = rep(50,5))
#' print(myga@solution)
mcga2 <- function(fitness, ...,
        min, max, population = gaControl("real-valued")$population,
        selection = gaControl("real-valued")$selection,
        crossover = byte_crossover, 
        mutation = byte_mutation,
        popSize = 50, 
        pcrossover = 0.8, 
        pmutation = 0.1, 
        elitism = base::max(1, round(popSize*0.05)), 
        maxiter = 100,
        run = maxiter,
        maxfitness = Inf,
        names = NULL,
        parallel = FALSE,
        monitor = gaMonitor,
        seed = NULL){


myga <- ga(type = "real-valued", fitness = fitness, ..., min = min, max = max, population = population, selection = selection, crossover = crossover, mutation = mutation, popSize = popSize, pcrossover = pcrossover, pmutation = pmutation, elitism = elitism, maxiter = maxiter, run = run, maxfitness = maxfitness, names = names, parallel = parallel, monitor = monitor, seed = seed)

return(myga)
}




