library(tsModel)

x <- rep(1:5, 2)
f <- gl(2, 5)

Lag(x, 1, f)
Lag(x, -1, f)
Lag(x, 0, f)

x <- rep(1:5, each = 2)
f <- gl(5, 2)

Lag(x, 1, f)
Lag(x, -1, f)
Lag(x, 0, f)

