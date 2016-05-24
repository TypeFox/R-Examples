
library("mlt")
set.seed(29)

n <- 20
### design has rank deficit; just for interface checking
### we need something better!
d <- data.frame(x1 = 1:n, x2 = 1:n + 1, y = rnorm(n))
m <- ctm(polynomial_basis(numeric_var("y", support = range(d$y)),
                            coef = c(TRUE, TRUE), ci = c(-Inf, 0)),
           shift = ~ x1 + x2, data = d)
mod <- mlt(m, data = d)

p <- predict(mod, newdata = d)

(p0 <- predict(mod$model$model, 
    newdata = expand.grid(d), coef = coef(mod)))
(p1 <- predict(mod, newdata = as.list(d)))
(p2 <- predict(mod, newdata = d, q = d$y[1]))

max(abs(p0 - as.vector(p1)))

all.equal(p1[cbind(1:n, 1:n, 1), drop = TRUE],
          drop(p2))

all.equal(p1[cbind(1:n, 1:n, 1:n), drop = TRUE],
          drop(p), check.attributes = FALSE)

predict(mod, newdata = list(x1 = 1:3, x2 = 2:3), p = c(.25, .5), type = "quantile")

simulate(mod, nsim = 1, seed = 291, interpolate = FALSE)

d$y <- gl(3, 1, ordered = TRUE)[rep(1:3, length = n)]

r <- as.basis(d$y) #as.basis(~ y, data = d, remove_intercept = TRUE,
#              contrasts.arg = list(y = function(n)
#                  contr.treatment(n, base = 3)),
#              ui = diff(diag(2)), ci = 0)

mod2 <- mlt(ctm(r, shift = ~ x1 + x2, data = d), data = d)

predict(mod2, q = unique(d$y))

predict(mod2, p = 1:9 / 10, type = "quantile")

simulate(mod2, nsim = 3, seed = 29)

predict(mod2, q = unique(d$y), type = "density")

predict(mod2, list(y = unique(d$y), x1 = 1:3, x2 = 2:3), type = "density")

