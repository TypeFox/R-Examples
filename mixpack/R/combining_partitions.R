#' Build a hierchical partition from posterior probabilities
#' 
#' This function applies the methodology described in [citar article]
#' to build a hierarchy of classes using the weights or probabilities 
#' that an element belongs to each class
#' @param post dataframe of probabilities/weights (\code{tau} must be strictly positive)
#' 
#' @param omega String giving the function name used to build the hierarchy. Available 
#' functions are: entr, prop, dich
#' 
#' @param lambda String giving the function name used to build the hierarchy. Available 
#' functions are: entr, demp, demp.mod, coda, coda.norm, prop
#' 
#' @param f_omega function with two parameters (\code{v_tau}, \code{a}). Parameter 
#' \code{v_tau} is a vector of probabilities, parameter \code{a} is the a selected class.
#' \code{omega}(\code{v_tau}, \code{a}) gives the representativeness of element with
#' probabities \code{v_tau} to class \code{a}
#' 
#' @param f_lambda function with three parameters (\code{v_tau}, \code{a}, \code{b}).
#' Parameter \code{v_tau} is a vector of probabilities, parameters \code{a} and \code{b}
#' are classes to be combined.
#' @export
get_hierarchical_partition <- function(post, omega, lambda, f_omega = NULL, f_lambda = NULL) {
  if( !missing(omega) & !missing(lambda) ){
    get_hierarchical_partition_cpp(post, omega, lambda)
  }else{
    if(is.null(f_omega) | is.null(f_lambda)){
      stop("Either omega and lambda or f_omega and f_lambda need to be defined")
    }
    if(!is.function(f_omega) | !is.function(f_lambda)){
      stop("Both f_omega and f_lambda needs to be a function")
    }
    if(length(formals(f_omega)) != 2){
      stop("f_omega needs two parameters. Probability vector and one position")
    }
    if(length(formals(f_lambda)) != 3){
      stop("f_lambda needs two parameters. Probability vector and two positions")
    }
    f_res = f_omega(c(0.5,0.5), 1)
    if(!is.numeric(f_res) & length(f_res) == 1){
      stop("f_omega needs to return a single number")
    }
    f_res = f_lambda(c(0.5,0.5), 1, 2)
    if(!is.numeric(f_res) & length(f_res) == 1){
      stop("f_lambda needs to return a single number")
    }
    get_hierarchical_partition_generic(post, f_omega, f_lambda)
  }
}

get_hierarchical_partition_cpp <- function(post, omega, lambda) {
  if(!omega %in% c('cnst', 'prop', 'dich')){
    stop(sprintf("Omega function %s is not available", omega))
  }
  if(!lambda %in% c('entr', 'demp', 'demp.mod', 'coda', 'coda.norm', 'prop')){
    stop(sprintf("Lambda function %s is not available", lambda))
  }
  if( (num <- sum(post == 0)) > 0){
    if( lambda %in% c('coda', 'coda.norm')){
      message(sprintf("%d zeros were replace by minimum machine number %e", num, .Machine$double.xmin))
      post[post==0] = .Machine$double.xmin
    }
  }
  .Call('mixpack_get_hierarchical_partition_cpp', PACKAGE = 'mixpack', post, omega, lambda)
}

get_hierarchical_partition_generic <- function(tau, omega, lambda) {
  ctau <- tau
  K <- ncol(ctau)
  partitions <- list()
  partitions[[K]] <- as.list(1:K)
  names(partitions[[K]]) <- sapply(partitions[[K]], part_name)
  for (k in K:2) {
    COMB <- t(expand.grid(1:k, 1:k))
    COMB <- COMB[, COMB[1, ] != COMB[2, ]]
    rownames(COMB) <- c("a", "b")
    colnames(COMB) <- col.names <- apply(COMB, 2, paste, collapse = "-")
    to_merge <- which.max(v <- apply(COMB, 2, function(ind) {
      a <- ind[1]
      b <- ind[2]
      sum(apply(ctau, 1, function(v_tau) omega(v_tau, a) * lambda(v_tau, a, b)))/sum(apply(ctau, 1, function(v_tau) omega(v_tau, 
        a)))
    }))
    part <- COMB[, to_merge]
    partitions[[k - 1]] <- b_absorbes_a(partitions[[k]], part["a"], part["b"])
    ctau[, part["b"]] <- ctau[, part["a"]] + ctau[, part["b"]]
    ctau <- ctau[, -part["a"]]
  }
  class(partitions) <- "hpartition"
  partitions
}

#' Merging components step
#' 
#' @param post Matrix with the posterior probabilities
#' @param omega String giving the function name used to build the hierarchy. Available 
#' functions are: entr, prop, dich
#' 
#' @param lambda String giving the function name used to build the hierarchy. Available 
#' functions are: entr, demp, demp.mod, coda, coda.norm, prop
#' 
#' @param f_omega function with two parameters (\code{v_tau}, \code{a}). Parameter 
#' \code{v_tau} is a vector of probabilities, parameter \code{a} is the a selected class.
#' \code{omega}(\code{v_tau}, \code{a}) gives the representativeness of element with
#' probabities \code{v_tau} to class \code{a}
#' 
#' @param f_lambda function with three parameters (\code{v_tau}, \code{a}, \code{b}).
#' Parameter \code{v_tau} is a vector of probabilities, parameters \code{a} and \code{b}
#' are classes to be combined.
#' @return partition returns a matrix with all values for all possible mergings using functions `omega` and `lambda`
#' @export
merge_step = function(post, omega, lambda, f_omega = NULL, f_lambda = NULL){
  if( !missing(omega) & !missing(lambda) ){
    merge_step_cpp(post, omega, lambda)
  }else{
    if(is.null(f_omega) | is.null(f_lambda)){
      stop("Either omega and lambda or f_omega and f_lambda need to be defined")
    }
    if(!is.function(f_omega) | !is.function(f_lambda)){
      stop("Both f_omega and f_lambda needs to be a function")
    }
    if(length(formals(f_omega)) != 2){
      stop("f_omega needs two parameters. Probability vector and one position")
    }
    if(length(formals(f_lambda)) != 3){
      stop("f_lambda needs two parameters. Probability vector and two positions")
    }
    f_res = f_omega(c(0.5,0.5), 1)
    if(!is.numeric(f_res) & length(f_res) == 1){
      stop("f_omega needs to return a single number")
    }
    f_res = f_lambda(c(0.5,0.5), 1, 2)
    if(!is.numeric(f_res) & length(f_res) == 1){
      stop("f_lambda needs to return a single number")
    }
    merge_step_generic(post, f_omega, f_lambda)
  }
}

merge_step_cpp <- function(post, omega, lambda) {
  if(!omega %in% c('cnst', 'prop', 'dich')){
    stop(sprintf("Omega function %s is not available", omega))
  }
  if(!lambda %in% c('entr', 'demp', 'demp.mod', 'coda', 'coda.norm', 'prop')){
    stop(sprintf("Lambda function %s is not available", lambda))
  }
  if( (num <- sum(post == 0)) > 0){
    if( lambda %in% c('coda', 'coda.norm')){
      message(sprintf("%d zeros were replace by minimum machine number %e", num, .Machine$double.xmin))
      post[post==0] = .Machine$double.xmin
    }
  }
  .Call('mixpack_merge_step_cpp', PACKAGE = 'mixpack', post, omega, lambda)
}

merge_step_generic <- function(tau, omega, lambda) {
  message("To be implenented")
}

#' Build a hierchical partition randomly from given K
#' 
#' This function return a hierachical partition contructed randonmly.
#' 
#' @param K number of initial groups
#' 
#' @export
get_random_hierarchical_partition <- function(K) {
  partitions <- list()
  partitions[[K]] <- as.list(1:K)
  names(partitions[[K]]) <- sapply(partitions[[K]], part_name)
  for (k in K:2) {
    COMB <- t(expand.grid(1:k, 1:k))
    COMB <- COMB[, COMB[1, ] != COMB[2, ]]
    rownames(COMB) <- c("a", "b")
    colnames(COMB) <- col.names <- apply(COMB, 2, paste, collapse = "-")
    to_merge <- sample(1:ncol(COMB), 1)
    part <- COMB[, to_merge]
    partitions[[k - 1]] <- b_absorbes_a(partitions[[k]], part["a"], part["b"])
  }
  class(partitions) <- "hpartition"
  partitions
}

#' Build a hierchical partition randomly from given K
#' 
#' This function return a hierachical partition contructed randonmly.
#' 
#' @param post posterior probability matrix
#' @param a first component to merge
#' @param b second component to merge
#' @return a matrix of posterior probabilities where components a and b are merged
#' 
#' @export
merge_components = function(post, a, b){
  if(!(1 <= a & a <= NCOL(post))){
    stop(sprintf('position a has to be lower or equal than %d', NCOL(post)))
  }
  if(!(1 <= b & b <= NCOL(post))){
    stop(sprintf('position b has to be lower or equal than %d', NCOL(post)))
  }
  A = a - 1
  B = b - 1
  .Call('mixpack_mergeComponents', PACKAGE = 'mixpack', post, A, B)
}

part_name <- function(part) sprintf("(%s)", paste(sort(part), collapse = ","))

## PART B absorbes A. In this function part A is incorporated to B and after that partition A is eliminated
b_absorbes_a <- function(partition, partA, partB) {
  if (partA == partB) {
    stop("Same part A and B")
  }
  if (!(partA %in% 1:length(partition) & partB %in% 1:length(partition))) {
    stop("Some part out of range")
  }
  new_partition <- partition
  new_partition[[partB]] <- c(new_partition[[partA]], new_partition[[partB]])
  new_partition[[partA]] <- NULL
  names(new_partition) <- sapply(new_partition, part_name)
  new_partition
}
#' Create a cluster from a partition
#' 
#' Given a matrix of tau and a partition decide in which part is classified each observation
#' @param tau matrix of posterioris
#' @param partition list of vectors containing the partition
#' @export
cluster_partition <- function(tau, partition) names(partition)[apply(do.call("cbind", lapply(partition, function(part) {
  if (is.vector(tau[, part])) 
    return(tau[, part])
  apply(tau[, part], 1, sum)
})), 1, which.max)] 

prop_partition_mult = function(tau, partition)
  do.call('cbind', lapply(partition, function(part){
    if(is.vector(tau[,part])) return(tau[,part])
    apply(tau[,part], 1, prod)^(1/length(part))
  }))
