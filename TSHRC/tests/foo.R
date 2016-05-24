
 library(TSHRC)

 n <- 200

 set.seed(42)

 x1 <- rexp(n)
 y1 <- rexp(n, rate = 1 / 2)
 t1 <- pmin(x1, y1)
 d1 <- as.numeric(t1 == x1)
 g1 <- rep(0, n)

 a <- 1.25
 b <- 1 / gamma(1 + 1 / a)
 x2 <- rweibull(n, shape = a, scale = b)
 y2 <- rexp(n, rate = 1 / 2)
 t2 <- pmin(x2, y2)
 d2 <- as.numeric(t2 == x2)
 g2 <- rep(1, n)

 twostage(c(t1, t2), c(d1, d2), c(g1, g2), nboot = 2500)

 twostage(c(t1, t2), c(d1, d2), c(g1, g2), nboot = 2500)
 
 twostage(c(t1, t2), c(d1, d2), c(g1, g2), nboot = 2500)

