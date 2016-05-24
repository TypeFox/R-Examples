
 library(pooh)

 set.seed(42)

 m <- 200
 n <- 500
 
 d <- seq(1:m)
 foo <- sample(d, n, replace = TRUE)
 bar <- sample(d, n, replace = TRUE)
 baz <- pmin(foo, bar)
 qux <- pmax(foo, bar)

 quacks <- tsort(baz, qux, d, strict = FALSE)
 try(tsort(baz, qux, d))
 any(baz == qux)

 length(quacks) == length(d)
 identical(sort(quacks), d)
 idx <- match(baz, quacks)
 jdx <- match(qux, quacks)
 all(idx <= jdx)

 quacks <- try(tsort(1:2, 2:1))

 tsort(character(0), character(0))
 tsort(character(0), character(0), domain = 1:5)

 try(tsort(1:5, 2:6))
 try(tsort(c(1:5, 1), c(2:6, 1)))
 try(tsort(c(1:5, 1), c(2:6, 1), strict = FALSE))

