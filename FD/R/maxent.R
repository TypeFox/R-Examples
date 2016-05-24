maxent <- function(constr, states, prior, tol = 1e-07, lambda = FALSE){
     
     # check input
    if (is.vector(constr)){
        means.names <- names(constr)
        constr <- matrix(constr, 1, length(constr) ) ; dimnames(constr) <- list("set1", means.names)
    }
    if (is.data.frame(constr)) constr <- as.matrix(constr)
    if (!is.numeric(constr)) stop("constr must only contain numeric values\n")
    if (!is.numeric(tol)) stop("tol must be numeric\n")
    if (!is.logical(lambda)) stop("lambda must be logical\n")
    if (is.vector(states)){
          s.names <- names(states)
          states <- matrix(states, nrow = 1, ncol = length(states) )
          dimnames(states) <- list("constraint", s.names)
         }
    if (is.data.frame(states)) states <- as.matrix(states)
    if (!is.numeric(states)) stop("states must only contain numeric values\n")
    if (any(states <= 0)) stop("states cannot contain zero or negative values\n")
    if (dim(states)[2] == 1 && dim(states)[1] > 1) states <- t(states)
    s.names <- dimnames(states)[[2]]
    c.names <- dimnames(states)[[1]]
    set.names <- dimnames(constr)[[1]]
    n.states <- dim(states)[2]
    n.traits <- dim(states)[1]
    n.sets <- dim(constr)[1]
    if (n.traits != dim(constr)[2]) stop("number of constraints in constr should be equal to number of constraints in states\n")
    if (missing(prior) ) {prior <- matrix(1 / n.states, n.sets, n.states)   ; dimnames(prior) <- list(set.names, s.names)}
    if (is.vector(prior)){
        if (length(prior) != n.states) {stop("number of states in prior should be equal to number in states\n")}
        if (n.sets == 1) {prior <- matrix(prior, 1, length(prior) ) ; dimnames(prior) <- list("set1", s.names)}
        else {prior <- matrix(rep(prior, n.sets), n.sets, length(prior), byrow = T ) ; dimnames(prior) <- list(set.names, s.names)}
    }     
    if (is.data.frame(prior)) prior <- as.matrix(prior)
    if (!is.numeric(prior)) stop("prior must only contain numeric values\n")
    if (dim(prior)[2] == 1 && dim(prior)[1] > 1) prior <- t(prior)
    if (dim(prior)[2] != n.states) stop("number of states in prior should be equal to number in states\n")
    if (dim(prior)[1] > 1 && dim(prior)[1] != n.sets) stop("number of rows in prior should be 1 or equal to number of rows in constr\n")
    if (dim(prior)[1] == 1) prior <- matrix(rep(prior[1,], n.sets), n.sets, ncol(prior), byrow = T )
    if (any(prior < 0 | prior > 1)) stop("prior must contain probabilities between 0 and 1\n")
    prior <- t(apply(prior, 1, function(x) x / sum(x))); dimnames(prior) <- list(set.names, s.names)
    if (any(is.na(constr)) || any(is.na(states)) || any(is.na(prior)) ) stop("no NA's allowed\n")
    
    # create storage objects
    allprobs <- matrix(NA, n.sets, n.states); dimnames(allprobs) <- list(set.names, s.names)
    moments <- matrix(NA, n.sets, n.traits) ; dimnames(moments) <- list(set.names, c.names)
    entropy <- rep(NA, n.sets) ; names(entropy) <- set.names
    iter <- rep(NA, n.sets) ; names(iter) <- set.names
    if (lambda){
       lambdas <- matrix(NA, n.sets, n.traits + 2)
       dimnames(lambdas) <- list(set.names, c("intercept", c.names, "prior") )
    }
    # FORTRAN loop
    for (i in 1:n.sets){
       itscale <- .Fortran("itscale5", as.double(t(states)), as.integer(n.states), as.integer(n.traits), as.double(constr[i,]), as.double(prior[i,]), prob = double(n.states), entropy = double(1), niter = integer(1), as.double(tol), moments = double(n.traits), PACKAGE = "FD")
       allprobs[i, ] <- itscale$prob
       moments[i, ] <- itscale$moments
       entropy[i] <- itscale$entropy
       iter[i] <- itscale$niter
       if (lambda) lambdas[i, ] <- coef(lm(log(itscale$prob) ~ t(states) + log(prior[i,]) ) )
   }

  # output
  res <- list()
  if (n.sets == 1){
     res$prob <- allprobs[1, ]
     res$moments <- moments[1, ]
     names(entropy) <- NULL
     res$entropy <- entropy
     names(iter) <- NULL
     res$iter <- iter
     if (lambda) res$lambda <- lambdas[1, ]
     res$constr <- constr[1, ]
     res$states <- states
     res$prior <- prior[1, ]
  }
  else{
    res$prob <- allprobs
    res$moments <- moments
    res$entropy <- entropy
    res$iter <- iter
    if (lambda) res$lambda <- lambdas
    res$constr <- constr
    res$states <- states
    res$prior <- prior
  }
  return(res)
}

