
library(modeltools)

d <- data.frame(x = rnorm(100), y = rnorm(100), z = runif(100))
d[["x"]][1:10] <- NA

a <- linearModel@dpp(y ~ x + z - 1, data = d, na.action = na.pass)
b <- na.omit(a)
mod1 <- linearModel@fit(b)

mod2 <- lm(y ~ x + z - 1, data = d)

nd <- data.frame(x = rnorm(100), z = runif(100))

stopifnot(identical(mod1$predict_response(nd), predict(mod2, newdata = nd)))

stopifnot(identical(coef(mod1), coef(mod2)))

u <- linearModel@fit
system.time(for (i in 1:100) mod1 <- u(b))
system.time(for (i in 1:100) mod2 <- lm(y ~ x + z - 1, data = d))

dn <- data.frame(x = rnorm(100), y = rnorm(100), z = runif(100))
all.equal(Predict(mod1, dn), Predict(mod2, dn))

system.time(for (i in 1:100) p1 <- Predict(mod1, dn))
system.time(for (i in 1:100) p2 <- Predict(mod2, dn))

system.time(for (i in 1:100) p1 <- predict(mod1, dn))
system.time(for (i in 1:100) p2 <- predict(mod2, dn))

### check bug fix: non-misssing `data' argument
df <- data.frame(y = 1:10, x = 1:10 + 1, z = 1:10 + 2)
mf <- ModelEnvFormula(y ~ x, data = df, other = list(part = ~ z))
stopifnot(isTRUE(all.equal(mf@get("part")$z, df[["z"]])))             
df2 <- df + 1
stopifnot(isTRUE(all.equal(mf@get("part", data = df2)$z, df2[["z"]])))

### ~ 1
df <- data.frame(y = 1:10)
mf <- ModelEnvFormula(y ~ 1, data = df)
x <- mf@get("designMatrix")
stopifnot(nrow(x) == 10 && all(x[,1] == 1))

### bugfix: subset was not correctly interpreted in `frame'
tmp <- function(formula, data = list(), subset = NULL) 
    ModelEnvFormula(formula, data, subset = subset, frame = parent.frame())
foo <- function(x, y, subset, ...) tmp(y ~ x, subset = subset, ...)
a <- 1:10     
b <- 1:10     
stopifnot(identical(foo(a, b, subset = 1:5)@get("response")[[1]],1:5)) 

x <- 1
y <- 2   
stopifnot(identical(foo(a, b, subset = 1:5)@get("response")[[1]],1:5))   

### subset problems
menv <- ModelEnvFormula(Species ~ ., data = iris, 
                        subset = (iris$Species != "virginica"))
stopifnot(nrow(menv@get("input")) == 100)
stopifnot(nrow(menv@get("input", data = iris)) == 150)

menv <- ModelEnvFormula(Species ~ ., data = iris, 
                        subset = (iris$Species != "virginica"), 
                        keep.subset = TRUE)
stopifnot(nrow(menv@get("input")) == 100)
stopifnot(nrow(menv@get("input", data = iris)) == 150)


###**********************************************************

stopifnot(!empty(menv))
menv1 <- new("ModelEnv")
stopifnot(empty(menv1))

### fixed in 0.2-17
dpp(linearModel, Sepal.Length ~ 1, data = iris, na.action = na.omit)
