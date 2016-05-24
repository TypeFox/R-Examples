#' Genetic Algorithm optimization
#' 
#' Function optimization through Genetic Algorithms.
#' 
#' Given a real-based or permutation-based function, and the associated search space,
#' \code{gaoptim} will perform a function maximization using the Genetic Algorithm
#' approach. For better performance, a real-number encoding is used.
#' 
#' All you need to get started is to provide a function and the associated search
#' space - there are sensible defaults to all the other parameters. On the other hand,
#' you can provide custom genetic-operators to control how your population will
#' \code{reproduce} and \code{mutate} (see the examples).
#' 
#' After setting the algorithm parameters, you can evolve your population and check the
#' results. You don't need to do this in one step, you can always evolve a small number
#' of generations and query the best solution found. If this solution doesn't fit your
#' needs, you can keep evolving your population - this approach saves time and computer
#' resources.
#' 
#' @name gaoptim
#' @docType package
#' @references Randy L. Haupt, Sue Ellen Haupt (2004). Practical genetic 
#' algorithms - 2nd ed.
#' @references Michalewicz, Zbigniew. Genetic Algorithms + Data Structures = Evolution
#' Programs - 3rd ed.
#' 
#' @references Luke, Sean. Department of Computer Science. George Mason University. 
#' Essentials of Metaheuristics - online version 1.2, July 2011. 

NULL

