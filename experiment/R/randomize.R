randomize <- function(data, group = c("Treat", "Control"), ratio =
                      NULL, indx = NULL, block = NULL, n.block = NULL,
                      match = NULL, complete = TRUE){  

  ## call
  call <- match.call()
  ## data 
  m <- length(group)
  if ((!is.null(call$block)) && (!is.null(call$match))) {
    stop("invalid inputs for `block' and `match'.")
  } else if (!is.null(call$block)) { ## blocking
    if ("formula" %in% class(block)) {
      tm <- terms(block)
      attr(tm, "intercept") <- 0
      data <- model.frame(tm, data = data, na.action = na.fail)
      X <- model.matrix(tm, data = data)
      block <- mahalanobis(X, apply(X, 2, mean), var(X))
    } else {
      block <- eval(call$block, envir = data)
    }
  } else { ## matching
    if (m != 2)
      stop("2 groups are required for matching.")
    if (!is.null(ratio))
      warning("`ratio' will be ignored.")
    if (complete) {
      if ("formula" %in% class(match)) {
        tm <- terms(match)
        attr(tm, "intercept") <- 0
        data <- model.frame(tm, data = data, na.action = na.fail)
        X <- model.matrix(tm, data = data)
        match <- mahalanobis(X, apply(X, 2, mean), var(X))
      } else {
        match <- eval(call$match, data)
      }
    } else {
      stop("`complete' should be TRUE for matching.")
    }
  }
    
  ## getting index
  n <- nrow(data)
  if (!is.null(indx))
    indx <- eval(indx, data)
  else if (is.null(rownames(data)))
    indx <- 1:n
  else
    indx <- rownames(data)
  
  ## groups
  if (is.null(ratio))
    ratio <- rep(1/m, m)
  ratio <- ratio/sum(ratio)
  if (sum(ratio < 0) > 1)
    stop("invalid input for `size'.")

  ## output
  res <- list(call = call, ratio = ratio)

  ## blocking and matching variable
  if (is.null(block) && is.null(match)) { 
    if (complete) { # complete randomization
      tmp <- ratio2size(n, ratio, group)
      ttt <- sample(tmp$vector, n, replace = FALSE)
    } else { # simple randomization
      ttt <- sample(group, n, replace = TRUE, prob = ratio)
    }
    names(ttt) <- indx
  } else if (is.null(match)) { ## blocking
    block.id <- rep(NA, n)
    if (is.null(n.block)) {
      tmp <- unique(block)
      n.block <- length(tmp)
      for (i in 1:n.block)
        block.id[block == tmp[i]] <- i 
    } else {
      tmp <- quantile(block, (0:(n.block-1))/n.block)
      for (i in 1:n.block)
        block.id[block >= tmp[i]] <- block.id[block >= tmp[i]] + 1
      if (sum(table(block.id) < m) > 0)
        stop("some blocks have too few observations.")
    }
    ttt <- rep(NA, n)
    names(ttt) <- names(block.id) <- indx
    for (i in 1:n.block) {
      howmany <- sum(block.id == i)
      if (complete) { # comlete randomization
        tmp <- ratio2size(howmany, ratio, group)
        ttt[block.id == i] <- sample(tmp$vector, howmany, replace = FALSE)
      } else {
        ttt[block.id == i] <- sample(group, sum(block.id == i),
            replace = TRUE, prob = ratio)
      }
    }
    res$block <- block
    res$block.id <- block.id
  } else { ## matching
    match.id <- ttt <- rep(NA, n)
    names(match.id) <- names(ttt) <- indx
    counter <- 1
    while (sum(is.na(match.id)) > 1) {
      unit <- sample(indx[is.na(match.id)], 1, replace = FALSE)
      diff <- abs(match[is.na(match.id)]-match[unit])
      mindiff <- names(sort(diff[diff>0]))[1]
      match.id[unit] <- match.id[mindiff] <- counter
      tmp <- sample(group, 2, replace = FALSE)
      ttt[unit] <- tmp[1]
      ttt[mindiff] <- tmp[2]
      counter <- counter + 1
    }
    res$match <- match
    res$match.id <- match.id
  }

  ## return the results
  res$treatment <- ttt
  res$data <- data
  class(res) <- "randomize"
  
  return(res)
}

###
### This converts ratio into size while randomly allocating remainders 
###

ratio2size <- function(n, ratio, group) {
  m <- length(ratio)

  size <- round(ratio * n, digits = 0)
  if (sum(size) > n) {
    tmp <- sample(1:length(size), sum(size)-n, replace = FALSE)
    size[tmp] <- size[tmp] - 1
  }
  if (sum(size) < n) {
    tmp <- sample(1:length(size), n-sum(size), replace = FALSE)
    size[tmp] <- size[tmp] + 1
  }
  allgroup <- NULL
  for (i in 1:m)
    allgroup <- c(allgroup, rep(group[i], size[i]))
  return(list(size = size, vector = allgroup))
  
}
