
 library(mcmc)

 x <- seq(0, 10, length = 10001)

 ### sub-exponentially light transformation

 b <- 0.5
 fsub <- morph(b = b)

 y <- unlist(Map(fsub$inverse, x))

 myfsub <- function(x) ifelse(x > 1 / b, exp(b * x) - exp(1) / 3,
     (x * b)^3 * exp(1) / 6 + x * b * exp(1) / 2)
 y2 <- myfsub(x)
 all.equal(y, y2, tolerance = 1e-14)

 z <- unlist(Map(fsub$transform, y))
 all.equal(z, x, tolerance = 1e-14)

 ### exponentially light transformation

 r <- 5 
 p <- 3
 fp3 <- morph(r = r)

 y <- unlist(Map(fp3$inverse, x))

 myfp3 <- function(x) ifelse(x < r, x, x + (x - r)^p)
 y2 <- myfp3(x)
 all.equal(y, y2, tolerance = 1e-14)

 z <- unlist(Map(fp3$transform, y))
 all.equal(z, x, tolerance = 1e-12)

 ### both together

 fboth <- morph(b = b, r = r)

 y <- unlist(Map(fboth$inverse, x))
 y2 <- myfsub(myfp3(x))
 all.equal(y, y2, tolerance = 1e-14)

 z <- unlist(Map(fboth$transform, y))
 all.equal(z, x, tolerance = 1e-12)

 ### exponentially light transformation with p != 3

 r <- 5 
 p <- 2.2
 fpo <- morph(r = r, p = p)

 y <- unlist(Map(fpo$inverse, x))

 myfpo <- function(x) ifelse(x < r, x, x + (x - r)^p)
 y2 <- myfpo(x)
 all.equal(y, y2, tolerance = 1e-14)

 z <- unlist(Map(fpo$transform, y))
 all.equal(z, x, tolerance = 1e-14)

