
library("basefun")

x <- 1:5
y <- as.double(1:4)
g <- gl(3, 1)
d <- expand.grid(x = x, y = y, g = g)
xvar <- numeric_var("x", support = x)

cb <- c(logx = log_basis(xvar, remove_intercept = TRUE), 
        X = as.basis(~ y + g, data = expand.grid(y = y, g = g)))

X <- model.matrix(cb, data = d)
stopifnot(nrow(X) == nrow(d))

p <- predict(cb, newdata = d, coef = rep(1, ncol(X)))
stopifnot(length(p) == nrow(d))

(p2 <- predict(cb, newdata = list(x = x, y = y, g = g), 
               coef = rep(1, ncol(X))))

stopifnot(all.equal(p, c(p2), check.attributes = FALSE))

(p4 <- predict(cb, newdata = mkgrid(cb, 4), coef = rep(1, ncol(X))))

stopifnot(all.equal(p, c(p4), check.attributes = FALSE))

p <- predict(cb, newdata = expand.grid(g = g, x = x, y = y), 
             coef = rep(1, ncol(X)))

(p2 <- predict(cb, newdata = list(x = x, y = y, g = g), 
               coef = rep(1, ncol(X)),
               dim = c(g = length(g), x = length(x), y = length(y))))
	
stopifnot(all.equal(p, c(p2), check.attributes = FALSE))

XX <- model.matrix(cb[["X"]], data = list(y = y, g = g), 
                   dim = c(y = length(y), g = length(g)))

pX <- predict(cb[["X"]], newdata = list(y = y, g = g), coef = rep(1, ncol(X) - 1))

pX2 <- predict(cb[["X"]], newdata = expand.grid(y = y, g = g), 
               coef = rep(1, ncol(X) - 1))

stopifnot(all.equal(c(pX), c(pX2), check.attributes = FALSE))

bb <- b(logx = log_basis(xvar, remove_intercept = TRUE),
        X = as.basis(~ y + g, data = expand.grid(y = y, g = g)))

X <- model.matrix(bb, data = d)
stopifnot(nrow(X) == nrow(d))  

p <- predict(bb, newdata = d, coef = rep(1, ncol(X)))
stopifnot(length(p) == nrow(d))

(p2 <- predict(bb, newdata = list(x = x, y = y, g = g), coef = rep(1, ncol(X))))

stopifnot(all.equal(p, c(p2), check.attributes = FALSE))

dd <- list(x = x, y = y[1:3], g = rep(g[1], 3))
(p3 <- predict(bb, newdata = dd, coef = rep(1, ncol(X)),
               dim = c(x = 5, y = 3, g = 1)))

stopifnot(all.equal(drop(p3), 
    matrix(p2[as.matrix(expand.grid(1:5, 1:3, 1))], nrow = 5), 
    check.attributes = FALSE))

