library(coop)
m <- 10
n <- 3

set.seed(1234)
x <- matrix(rnorm(m*n), m, n)
x[c(1,3,9), c(1, 1, 2)] <- NA

use <- "everything"
t1 <- cor(x, use=use)
t2 <- pcor(x, use=use)
stopifnot(all.equal(t1, t2))

use <- "all"
t1 <- try(cor(x, use=use), silent=TRUE)
t2 <- try(pcor(x, use=use), silent=TRUE)
stopifnot(inherits(t1, "try-error") && inherits(t2, "try-error"))

use <- "complete"
t1 <- cor(x, use=use)
t2 <- pcor(x, use=use)
stopifnot(all.equal(t1, t2))




set.seed(1234)
x <- rnorm(10)
y <- runif(10)
x[c(1,3,9)] <- NA

use <- "everything"
t1 <- cor(x, y, use=use)
t2 <- pcor(x, y, use=use)
stopifnot(all.equal(t1, t2))

use <- "all"
t1 <- try(cor(x, y, use=use), silent=TRUE)
t2 <- try(pcor(x, y, use=use), silent=TRUE)
stopifnot(inherits(t1, "try-error") && inherits(t2, "try-error"))

use <- "complete"
t1 <- cor(x, y, use=use)
t2 <- pcor(x, y, use=use)
stopifnot(all.equal(t1, t2))
