
 library(rcdd)

 d <- 5
 npoint <- 50
 set.seed(42)
 x <- matrix(rnorm(d * npoint), ncol = d)

 foo <- cbind(1, x)
 foo <- cbind(0, foo)
 foo <- round(1e4 * foo)
 bar <- rep(1e4, length(foo))
 baz <- qdq(foo, bar)

 baz[1:3, ]
 
 out <- scdd(baz, inputincidence = TRUE, representation = "V")
 names(out)

 length(out$inputincidence)
 inies <- sapply(out$inputincidence, length) > 0
 length(inies)
 sum(inies)

