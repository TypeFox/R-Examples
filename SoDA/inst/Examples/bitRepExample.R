m <- 24; x <- .54321
pwrs <- 2^(1:m)
xpwrs <- x*pwrs
xrep <- trunc(xpwrs)/pwrs
bits <- c(xrep[[1]], diff(xrep))*pwrs
bits
