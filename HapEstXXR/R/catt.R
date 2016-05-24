catt <-
  function(y, x, score = c(0, 1, 2)) {
    miss <- unique(c(which(is.na(y)), which(is.na(x))))
    n.miss <- length(miss)
    if(n.miss > 0) {
      y <- y[-miss]
      x <- x[-miss]
    }
    if(!all((y == 0) | (y == 1))) 
      stop("y should be only 0 or 1.")
    if(!all((x == 0) | (x == 1) |(x == 2))) 
      stop("x should be only 0, 1 or 2.")
    ca <- x [y == 1]
    co <- x [y == 0]
    htca <- table(ca)
    htco <- table(co)
    A <- matrix(0, 2, 3)
    colnames(A) <- c(0, 1, 2)
    rownames(A) <- c(0, 1)
    A[1, names(htca)] <- htca
    A[2, names(htco)] <- htco
    ptt <- prop.trend.test(A[1, ], colSums(A), score = score)
    res <- list("2x3-table" = A, 
                chisq = as.numeric(ptt$statistic), 
                df = as.numeric(ptt$parameter), 
                p.value = as.numeric(ptt$p.value), 
                n.miss = n.miss)
    return(res)
  }
