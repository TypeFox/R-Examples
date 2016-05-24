## convenience function to set up sampling function
## based on parameter specifications
make_rnorm_fun <- function(par) {
  m <- par[1L]
  s <- par[2L]
  function(n) rnorm(n, mean = m, sd = s)
}

make_sample_fun <- function(par) {
  x <- if(is.null(dim(par))) par[1L] else par[, 1L]
  p <- if(is.null(dim(par))) par[2L] else par[, 2L]
  if(isTRUE(all.equal(sum(p), 1))) {
    function(n) x[sample.int(length(x), size = n, prob = p, replace = TRUE)]
  } else {
    function(n) {
      np <- n * p / sum(p)
      stopifnot(isTRUE(all.equal(np, round(np))))
      rep(x, round(np))[sample.int(round(sum(np)))]
    }
  }
}



simRaschmix <- function(design, extremes = FALSE, attributes = TRUE,...) {

  if(!(is.list(design) && length(design) == 2L & all(names(design) %in% 
                               c("ability", "difficulty")))){
    ## character specification
    if(is.character(design)) {
      nobs <- 1800
      switch(match.arg(design, c("rost1", "rost2", "rost3")),
             ## Rost gives item easiness parameter which are converted to
             ## item difficulty parameters here.
             "rost1" = {
               weights <- 1
               ability <- array(c(cbind(c(-2.7, -0.9, 0.9, 2.7), rep(1, 4))),
                         dim = c(4, 2, 1))
               difficulty <- matrix(c(2.7,2.1,1.5,0.9,0.3,-0.3,-0.9,-1.5,-2.1,-2.7), ncol = 1)
             },
             "rost2" = {
               weights <- c(1,1)
               ability <- array(c(cbind(c(-2.7, -0.9, 0.9, 2.7), rep(1, 4)),
                                  cbind(c(-2.7, -0.9, 0.9, 2.7), rep(1, 4))),
                                dim = c(4, 2, 2))
               difficulty <- cbind(c(2.7,2.1,1.5,0.9,0.3,-0.3,-0.9,-1.5,-2.1,-2.7),
                                   c(-2.7,-2.1,-1.5,-0.9,-0.3,0.3,0.9,1.5,2.1,2.7))
             },
             "rost3" = {
               weights <- c(4,2,3)
               ab.c3 <- numeric(4)
               ab.c3[sample.int(4, size = 1, prob = rep(1,4)/4)] <- 1
               ability <- array(c(cbind(c(-2.7, -0.9, 0.9, 2.7), rep(1, 4)),
                                  cbind(c(-2.7, -0.9, 0.9, 2.7), rep(1, 4)),
                                  cbind(c(-2.7, -0.9, 0.9, 2.7), ab.c3)),
                                dim = c(4, 2, 3))
               difficulty <- cbind(c(2.7,2.1,1.5,0.9,0.3,-0.3,-0.9,-1.5,-2.1,-2.7),
                                   c(-2.7,-2.1,-1.5,-0.9,-0.3,0.3,0.9,1.5,2.1,2.7),
                                   c(-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5))
             })
    } else if(is.list(design) & length(design) == 4L & all(names(design) %in% c("ability", "difficulty", "nobs", "weights"))) {
      nobs <- design$nobs
      weights <- design$weights
      ability <- design$ability
      difficulty <- design$difficulty
    }
    ## process weights, can be in three formats:
    ## (1) function(n)
    ## (2) vector of probabilities (summing to 1)
    ## (3) vector of integer weights (summing to n, or an integer division thereof)

    ## reduce (2)/(3) to case (1)
    if(!is.function(weights)) weights <- make_sample_fun(cbind(seq_along(weights), weights))
    ## set up weights
    cluster <- weights(nobs)
    weights <- table(cluster)
    k <- length(weights)
    
    ## ability specification (in three clusters?), can be in three formats:
    ## (1) list(k) of function(n)
    ## (2) matrix(2, k) with means and standard deviations per cluster
    ## (3) array(., 2, k) with abilities per weights/probabilities per cluster
    
    ## reduce (2)/(3) to case (1)
    if(!(is.list(ability) && all(sapply(ability, is.function)))) {
      if(length(dim(ability)) == 2L & dim(ability)[1L] == 2L & dim(ability)[2L] == k) {
        ability <- lapply(1:k, function(i) make_rnorm_fun(ability[, i]))
      } else if(length(dim(ability)) == 3L & dim(ability)[2L] == 2L & dim(ability)[3L] == k) {
        ability <- lapply(1:k, function(i) make_sample_fun(ability[, , i]))  
      } else {
        stop("unknown ability specification")
      }
    }
    ## set up ability within each cluster
    ab <- numeric(nobs)
    for(i in 1:k) ab[cluster == i] <- ability[[i]](weights[i])
    ability <- ab
    
    ## difficulty specification, can be in 2 formats
    ## (1) matrix(n, .) with difficulties per person
    ## (2) matrix(., k) with difficulties per cluster
    
    if(!(length(dim(difficulty)) == 2L & dim(difficulty)[1L] == nobs)){
      if(length(dim(difficulty)) == 2L && dim(difficulty)[2L] == k){
        ## transform from per cluster to per person
        difficulty.o <- difficulty
        dif <- matrix(NA, nrow = nobs, ncol = dim(difficulty)[1L])
        for (i in 1:k) dif[cluster == i,] <- rep(difficulty[,i], each = weights[i])
        difficulty <- dif
      }
      else {
        stop("unknown difficulty specification")
      }
    }
  } else {
    ability <- design$ability
    difficulty <- design$difficulty
    cluster <- NULL
  }
  
  ## generate response
  stopifnot(NROW(ability) == NROW(difficulty))
  linpred <- ability - difficulty
  prob <- plogis(linpred)
  rval <- rbinom(n = length(prob), prob = prob, size = 1)
  dim(rval) <- dim(prob)
  
  ## remove observations with extreme raw scores
  if (!extremes){
    w <- which(rowSums(rval) == 0 | rowSums(rval) == ncol(rval))
    if (length(w) > 0){
      rval <- rval[-w,]
      ability <- ability[-w]
      cluster <- cluster[-w]
    }
  }
  
  ## attach attributes
  if(attributes) {
    attr(rval, "ability") <- ability
    #if(exists("difficulty.o")) attr(rval, "difficulty") <- difficulty.o
    attr(rval, "difficulty") <- if(exists("difficulty.o")) difficulty.o else difficulty
    attr(rval, "cluster") <- cluster
  }
  
  ## return
  return(rval)
}
