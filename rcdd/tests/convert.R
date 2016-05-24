
 library(rcdd)

 set.seed(42)
 x <- matrix(runif(4), 2, 2)
 print(x)

 y <- d2q(x)
 print(y)

 xfoo <- x
 attr(xfoo, "blurfle") <- "woof!"
 class(xfoo) <- c("foo", "bar", "baz")
 yfoo <- d2q(xfoo)
 print(yfoo)

 z <- q2d(y)
 print(z)

 x - z

 foo <- d2q(0.5)
 print(foo)
 bar <- q2q(foo)
 print(bar)
 bar[1] <- "-15/3"
 baz <- q2q(bar)
 print(baz)

 numer <- as.character(seq(-3, 3))
 denom <- as.character(seq(along = numer))
 qux <- z2q(numer, denom)
 print(qux)

 numer <- seq(-3, 3)
 denom <- seq(along = numer)
 qux <- z2q(numer, denom)
 print(qux)

 qux <- z2q(numer, denom, canonicalize = FALSE)
 print(qux)

