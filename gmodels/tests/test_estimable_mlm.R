library(gmodels)

y <- cbind(rnorm(100), rnorm(100))

x1 <- rnorm(100)
x2 <- rnorm(100)

cm <- t(matrix(c(0, 1,-1)))
lm.1 <- lm(y ~ x1 + x2)

estimable(lm.1, cm)

## >> Error in coef(object) : object 'object' not found

gmodels:::estimable.mlm(lm.1, cm)
