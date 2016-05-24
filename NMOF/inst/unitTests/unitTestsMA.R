## -*- truncate-lines: t; -*-
test.MA <- function() {
    x <- rnorm(100L); myMA <- numeric(length(x)); order <- 5L
    for (i in order:length(x))
        myMA[i] <- sum(x[(i - order + 1):i])/order

    checkEquals(MA(x, order = order)[-(1:(order-1))],
                myMA[-(1:(order-1))])

    x <- rnorm(100L); myMA <- numeric(length(x)); order <- 1L
    checkEquals(x, MA(x, order = order))
    checkEquals(x, MA(x, order = order, pad = NA))
}
