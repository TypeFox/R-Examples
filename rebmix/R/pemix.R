pemix <- function(x = NULL, lower.tail = TRUE, log.p = FALSE, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  if (is.null(x)) {
    stop(sQuote("x"), " must not be NULL!", call. = FALSE)
  }

  if ((!is.numeric(x)) && (!is.data.frame(x))) {
    stop(sQuote("x"), " numeric or data frame is requested!", call. = FALSE)
  }

  x <- as.matrix(x)

  d <- ncol(x)  
  n <- nrow(x)
  
  y <- array(data = 0.0, dim = n, dimnames = NULL)

  if (lower.tail == TRUE) {
    for (i in 1:n) {
      y[i] <- sum(apply(x <= x[i,], 1, all))
    }
  }
  else {
    for (i in 1:n) {
      y[i] <- sum(apply(x > x[i,], 1, all))
    }
  }

  y <- y / n

  if (log.p == TRUE) {
    y <- log(y)
  }
  
  output <- as.data.frame(cbind(x, y), stringsAsFactors = FALSE)

  colnames(output) <- c(paste("x", if (d > 1) 1:d else "", sep = ""), "F")  

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  invisible(output)
} ## pemix
