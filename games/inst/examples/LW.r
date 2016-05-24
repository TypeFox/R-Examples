x <- rexp(10)
w <- LW(x)
all.equal(x, w * exp(w))
