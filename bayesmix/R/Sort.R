## Function taken from e1071
permutations <- function (n) {
    if (n == 1) 
        return(matrix(1))
    else if (n < 2) 
        stop("n must be a positive integer")
    z <- matrix(1)
    for (i in seq_len(n)[-1]) {
        x <- cbind(z, i)
        a <- c(seq_len(i), seq_len(i)[-i])
        z <- matrix(0, ncol = ncol(x), nrow = i * nrow(x))
        z[seq_len(nrow(x)), ] <- x
        for (j in seq_len(i)[-1] - 1) {
            z[j * nrow(x) + seq_len(nrow(x)), ] <- x[, a[seq_len(i) + j]]
        }
    }
    dimnames(z) <- NULL
    z
}

Sort <- function(x, by = NULL) {
  if (!(inherits(x, "jags") && inherits(x$model, "BMMmodel"))) 
    stop("Use only with 'jags' objects with model of class 'BMMmodel'.")
  x.old <- x
  n <- dim(x$results)
  if (is.null(by)) by <- x$variables
  else by <- x$variables[pmatch(by, x$variables)]
  by <- by[1]
  if (is.na(by)) stop("by not specified correctly")
  index <- grep(by, colnames(x$results))
  nn <- length(index)
  if (nn != x$model$data$k) stop("by not specified correctly")
  dd <- order(row(x$results[,index]), x$results[,index])
  ind <- apply(x$results[, index], 1, order)
  for (name in x$variables) {
    ii <- grep(name, colnames(x$results))
    if (length(ii) == nn) {
      x$results[,ii] <- matrix(x$results[,ii][dd], nrow = n[1], byrow = TRUE)
    }
    else if (length(levels(as.factor(x$results[,ii]))) == x$model$data$k) {
      ps <- permutations(x$model$data$k)
      for (j in seq_len(nrow(ps))) {
        ps1 <- ps[j,]
        index <- apply(ind, 2, function(x) all(x == ps1))
        if (any(index)) {
          dummy <- factor(x$results[index,ii], levels = seq_len(x$model$data$k))
          levels(dummy) <- order(ps1)
          x$results[index,ii] <- as.numeric(levels(dummy))[as.integer(dummy)]
        }
      }
    }
    else if (length(ii) != 1) {
      warning("Sorting not successful. Original object returned!")
      return(x.old)
    }
  }
  x
}

