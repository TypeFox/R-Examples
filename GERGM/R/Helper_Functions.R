# Directional derivative statistics next These are needed for the Gibbs sampler

# din2star
din2star <- function(i, j, net) {
  nodes <- nrow(net)
  others <- (1:nodes)[-c(i, j)]
  sum(net[cbind(others, j)])
}

# dout2star
dout2star <- function(i, j, net) {
  nodes <- nrow(net)
  others <- (1:nodes)[-c(i, j)]
  sum(net[cbind(i, others)])
}

# dedgeweight
dedgeweight = function(i, j) {
  1
}

# drecip
drecip <- function(i, j, net) {
  net[j, i]
}

# dctriads
dctriads <- function(i, j, net) {
  nodes <- nrow(net)
  others <- (1:nodes)[-c(i, j)]
  triples <- cbind(i, j, others)
  sum(net[triples[, c(2, 3)]] * net[triples[, c(3, 1)]])
}

# dttriads
dttriads <- function(i, j, net) {
  nodes <- nrow(net)
  others <- (1:nodes)[-c(i, j)]
  triples <- cbind(i, j, others)
  t2 <- sum(net[triples[, c(2, 3)]] * net[triples[, c(1, 3)]])
  t3 <- sum(net[triples[, c(3, 2)]] * net[triples[, c(3, 1)]])
  t4 <- sum(net[triples[, c(3, 2)]] * net[triples[, c(1, 3)]])
  return(t2 + t3 + t4)
}

# dtriads (undirected)
dtriads <- function(i, j, net){
  nodes <- nrow(net)
  others <- (1:nodes)[-c(i, j)]
  triples <- cbind(i, j, others)
  stat <- sum(net[triples[, c(2, 3)]] * net[triples[, c(1, 3)]])
  return(stat)
}


# dh function weight w_{i,j} will be conditioned upon
# Calculate the marginal change in the network
dh <- function(net, statistics, i, j) {
  temp <- c(dout2star(i, j, net), din2star(i, j, net), dctriads(i, j, net),
            drecip(i, j, net), dttriads(i, j, net), dedgeweight(i, j))
  if(length(temp) != length(statistics)){
    stop("Development ERROR! Please email mdenny@psu.edu! The dh() internal function in Helper_Functions.R has been supplied an incorrect number of statistics.")
  }
  value <- temp[statistics > 0]
  return(value)
}





# Functions needed for likelihood function calculations

dgam01 <- function(x, exbz, exlam) {
  x <- x + 0.1
  return(dgamma(x, shape = exbz / exlam, scale = exlam) /
           (1 - pgamma(0.1, shape = exbz / exlam, scale = exlam)))
}

pgam01 <- function(x, exbz, exlam) {
  x <- x + 0.1
  return((pgamma(x, shape = exbz / exlam, scale = exlam) -
            pgamma(0.1, shape = exbz / exlam, scale = exlam)) /
           (1 - pgamma(0.1, shape = exbz / exlam, scale = exlam)))
}

qgam01 <- function(p, exbz, exlam) {
  x <- qgamma(pgamma(0.1, shape = exbz / exlam, scale = exlam) +
                p * (1 - pgamma(0.1, shape = exbz / exlam, scale = exlam)),
              shape = exbz / exlam, scale = exlam)
  return(x - 0.1)
}

dst <- function(x, m, sig, df) {
  return(1 / sig * dt((x - m) / sig, df))
}

pst <- function(x, m, sig, df) {
  return(pt((x - m) / sig, df))
}

qst <- function(u, m, sig, df) {
  return(qt(u, df) * sig + m)
}

indeg <- function(net) {
  id <- numeric(nrow(net))
  for (i in 1:length(id)) {
    id[i] <- sum(net[, i])
  }
  return(id)
}

outdeg <- function(net) {
  od <- numeric(nrow(net))
  for (i in 1:length(od)) {
    od[i] <- sum(net[i, ])
  }
  return(od)
}

# Convert an observed network to edge weight vectors x and y

net2xy <- function(net, statistics, directed) {
  y <- NULL
  x <- NULL
  nodes <- nrow(net)
  if (directed == TRUE) {
    for (i in 1:nodes) {
      for (j in (1:nodes)[-i]) {
        y <- c(y, net[i, j])
        x <- rbind(x, dh(net, statistics, i, j))
      }
    }
  }
  if (directed == FALSE) {
    for (i in 1:nodes) {
      for (j in (1:nodes)[-i]) {
        y <- c(y, net[i, j])
        x <- rbind(x, dh(net, statistics, i, j))
      }
    }
#     for (i in 2:nodes) {
#       for (j in (1:(i - 1))) {
#         y <- c(y, net[i, j])
#         x <- rbind(x, dh(net, statistics, i, j))
#       }
#     }
  }
  return(list(y = y, x = x))
}

# Draw a random value either uniform or according to density
rtexp <- function(n, lambda) {
  # lambda is the scalar-valued parameter n is the number of draws
  u <- runif(n)
  if (lambda != 0) {
    temp = exp(lambda)
    temp[is.infinite(temp)] = 1e+32
    x <- log(1 + u * (temp - 1)) / lambda
    x
  } else u
}

# The conditional density of each weight from a sample
dtexp <- function(x, lambda) {
  den <- numeric(length(x))
  den[which(lambda != 0)] <- exp(x[which(lambda != 0)] *
                                   lambda[which(lambda != 0)]) /
    (1 / lambda[which(lambda != 0)] * (exp(lambda[which(lambda != 0)]) - 1))
  den[which(lambda == 0)] <- 1
  return(den)
}

# Function to generate dispersed unit interval
rdisp <- function(n) {
  pnorm(abs(rnorm(n) * 6))
}

# Log likelihood function calculations

# pseudolikelihood given theta#
pl <- function(theta, y, x) {
  return(sum(log(dtexp(y, x %*% theta))))
}

# add to console output field in GERGM_Object
#GERGM_Object <- store_console_output(GERGM_Object,addition)
store_console_output <- function(GERGM_Object,addition){
  if("list" %in% class(addition)){
    addition <- as.character(unlist(addition))
  }
  if(is.null(GERGM_Object@console_output)){
    GERGM_Object@console_output <- addition
  }else{
    GERGM_Object@console_output <- c(GERGM_Object@console_output,addition)
  }
  return(GERGM_Object)
}

# Parse a formula entry
parse_formula_term <- function(term,
                               possible_structural_terms,
                               possible_covariate_terms,
                               possible_network_terms){
  # split up the formula term
  parsed <- stringr::str_split(term,"[\\(\\)]+")[[1]]
  # get rid of trailing spaces
  rm_spaces <- which(parsed == "")
  if(length(rm_spaces) > 0){
    parsed <-  parsed[-rm_spaces]
  }

  # generate return list object
  return_list <- list(term = parsed[1],
                      weight = 1,
                      covariate = NA,
                      base = NA,
                      network = NA,
                      threshold = 0,
                      levels = NA,
                      same = NA,
                      parens_no_arg = NA,
                      network_matrix_object = NA,
                      num_levels = NA,
                      base_index = NA)
  possible_fields <- c("term","alpha","covariate", "base", "network",
                       "threshold", "levels", "same", "")
  # if there is an argument to the term -- this will be a lazy implementation
  # where if you do not get the name right, it will simply not be set and a
  # warning will be thrown.
  if(length(parsed) > 1){

    # take the arguments and further deparse them and turn them into a list
    # object, removing any quotes and further splitting each argument on an
    # equals sign so we can get the key value pair.
    args <- stringr::str_split(parsed[2],",")[[1]]
    args <- as.list(args)
    for(i in 1:length(args)){
      args[[i]] <- stringr::str_split(args[[i]],"=")[[1]]
      for(j in 1:length(args[[i]])){
        args[[i]][j] <- stringr::str_replace_all(args[[i]][j],"\"","")[[1]]
        args[[i]][j] <- stringr::str_replace_all(args[[i]][j],"\'","")[[1]]
      }
    }

    # if we are only supplied one term, without a name, we will default to
    # setting the exponential weight if it is a structural term and the
    # covariate name if it is a covariate term.

    ###### for structural terms ######
    if(return_list$term %in% possible_structural_terms){
      if(length(args[[1]]) == 1){
        # if an argument is supplied without an expicit argument name
        # assignment, then treat it as a weight and check if it is numeric
        # make sure the argument is numeric
        args[[1]] <- as.numeric(args[[1]])
        if(is.numeric(args[[1]])){
          return_list$weight <- args[[1]]
        }else{
          stop(paste("You must supply a numeric weight for structural covariate:", return_list$term))
        }
      }else{
        which_arg <- which(possible_fields == args[[1]][1])
        if(length(which_arg) == 0){
          stop(paste("You supplied an argument:",args[[1]][1],"which is not recognized."))
        }else{
          return_list[[which_arg]] <- args[[1]][2]
        }
      }
      ###### for covariates ######
    }else if(return_list$term %in% possible_covariate_terms){
      if(length(args[[1]]) == 1){
        # if an argument is supplied without an expicit argument name
        # assignment, then treat it as a weight and check if it is numeric
        if(is.character(args[[1]])){
          return_list$covariate <- args[[1]]
        }else{
          stop(paste("You must supply a valid name for each node level covariate. You specified:", return_list$term))
        }
      }else{
        which_arg <- which(possible_fields == args[[1]][1])
        if(length(which_arg) == 0){
          stop(paste("You supplied an argument:",args[[1]][2],"which is not recognized."))
        }else{
          return_list[[which_arg]] <- args[[1]][2]
        }
      }
      ###### for network covariates ######
    }else if(return_list$term %in% possible_network_terms){
      if(length(args[[1]]) == 1){
        # if an argument is supplied without an expicit argument name
        # assignment, then treat it as a weight and check if it is numeric
        if(is.character(args[[1]])){
          return_list$network <- args[[1]]
        }else{
          stop(paste("You must supply a valid name for each network covariate. You specified:", return_list$term))
        }
      }else{
        which_arg <- which(possible_fields == args[[1]][1])
        if(length(which_arg) == 0){
          stop(paste("You supplied an argument:",args[[1]][2],"which is not recognized."))
        }else{
          return_list[[which_arg]] <- args[[1]][2]
        }
      }
      # check to make sure a valid network matrix is specified
      if(!is.na(return_list$network) & !is.null(return_list$network)){
        temp_net <- dynGet(as.character(as.character(return_list$network)),
               ifnotfound = get(as.character(as.character(return_list$network))))
          return_list$network_matrix_object <- temp_net
        if(class(return_list$network_matrix_object)!= "matrix"){
          stop(paste("You must supply network covariates as matrix objects."))
        }
      }else{
        stop(paste("The network covariate matrix:",return_list$network,"does not appear to exist, please check that it is loaded in your current R session."))
      }
    }else{
      # throw an error because the term was not valid
      stop(paste("You specified term:",return_list$term,"which is not recognized. Please respecify."))
    }

    # if we are given two arguments.
    if(length(args) > 1){
      # loop through the rest of the arguments
      for(i in 2:length(args)){
        ###### for structural terms ######
        if(return_list$term %in% possible_structural_terms){
          which_arg <- which(possible_fields == args[[i]][1])
          if(length(which_arg) == 0){
            stop(paste("You supplied an argument:",args[[i]][2],"which is not recognized."))
          }else{
            return_list[[which_arg]] <- args[[i]][2]
          }
          ###### for covariates ######
        }else{
          which_arg <- which(possible_fields == args[[i]][1])
          if(length(which_arg) == 0){
            stop(paste("You supplied an argument:",args[[i]][2],"which is not recognized."))
          }else{
            return_list[[which_arg]] <- args[[i]][2]
          }
        }
      }
    }
  }

  # return everything
  return(return_list)
}
