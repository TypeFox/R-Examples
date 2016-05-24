
set.seed(290875)
library("party")

### get rid of the NAMESPACE
attach(asNamespace("party"))

### 
###
###    Regression tests for independence tests
###    
###    functions defined in file `./src/IndependenceTest.c'    

### tests for function C_IndependenceTest
xval <- rnorm(60)
x <- matrix(rank(xval), ncol = 1)
yf <- gl(2, 30)
y <- sapply(levels(yf), function(l) as.numeric(yf == l))[,1,drop = FALSE]
weights <- rep(1, nrow(x))

varctrl <- new("VariableControl")
lec <- new("LinStatExpectCovar", ncol(x), ncol(y))

res <- .Call("R_IndependenceTest", x, y, weights,  
              lec, varctrl, PACKAGE = "party")
print(res)
wmw <- wilcox.test(xval ~ yf, exact = FALSE, correct = FALSE)
print(wmw)
stopifnot(isequal(res[2], wmw$p.value))

xval <- rnorm(60)
x <- matrix(rank(xval), ncol = 1)
yf <- gl(3, 20)
y <- sapply(levels(yf), function(l) as.numeric(yf == l))
weights <- rep(1, nrow(x))

varctrl <- new("VariableControl")
varctrl@teststat <- factor("quad", levels = c("max", "quad"))
print(varctrl)
lec <- new("LinStatExpectCovarMPinv", ncol(x), ncol(y))
lec@rank <- 0
lec@MPinv <- matrix(0, nrow = ncol(x) * ncol(y), ncol = ncol(x) * ncol(y))
lec@svdmem <- new("svd_mem", ncol(x) * ncol(y))

res <- .Call("R_IndependenceTest", x, y, weights,  
              lec, varctrl, PACKAGE = "party")
print(res)
kw <- kruskal.test(xval ~ yf)
print(kw)
stopifnot(isequal(res[2], kw$p.value))
stopifnot(isequal(lec@rank, kw$parameter))

tmp <- x
x <- y
y <- tmp
varctrl <- new("VariableControl")
varctrl@teststat <- factor("quad", levels = c("max", "quad"))
print(varctrl)
lec <- new("LinStatExpectCovarMPinv", ncol(x), ncol(y))
lec@rank <- 0
lec@MPinv <- matrix(0, nrow = ncol(x) * ncol(y), ncol = ncol(x) * ncol(y))
lec@svdmem <- new("svd_mem", ncol(x) * ncol(y))

res <- .Call("R_IndependenceTest", x, y, weights,  
              lec, varctrl, PACKAGE = "party")
print(res)
kw <- kruskal.test(xval ~ yf)
print(kw)
stopifnot(isequal(res[2], kw$p.value))

xf <- gl(3, 2000)
x <- sapply(levels(xf), function(l) as.numeric(xf == l))
yf <- gl(3, 2000)[sample(1:6000)]
y <- sapply(levels(yf), function(l) as.numeric(yf == l))
weights <- rep(1, nrow(x))

varctrl <- new("VariableControl")
varctrl@teststat <- factor("quad", levels = c("max", "quad"))
print(varctrl)
lec <- new("LinStatExpectCovarMPinv", ncol(x), ncol(y))
lec@rank <- 0
lec@MPinv <- matrix(0, nrow = ncol(x) * ncol(y), ncol = ncol(x) * ncol(y))
lec@svdmem <- new("svd_mem", ncol(x) * ncol(y))

res <- .Call("R_IndependenceTest", x, y, weights,  
              lec, varctrl, PACKAGE = "party")
print(res)
chis <- chisq.test(table(xf, yf), correct = FALSE)
print(chis)
stopifnot(isequal(round(res[2],3), round(chis$p.value,3)))

### unbalanced data
xval <- rnorm(60)
x <- matrix(rank(xval), ncol = 1)
yf <- factor(rnorm(60) > 1)
y <- sapply(levels(yf), function(l) as.numeric(yf == l)) #[,1,drop = FALSE]
weights <- rep(1, nrow(x))

varctrl <- new("VariableControl")
lec <- new("LinStatExpectCovar", ncol(x), ncol(y))

res <- .Call("R_IndependenceTest", x, y, weights,
              lec, varctrl, PACKAGE = "party")
print(res)
wmw <- wilcox.test(xval ~ yf, exact = FALSE, correct = FALSE)
print(wmw)
stopifnot(isequal(res[2], wmw$p.value))

varctrl <- new("VariableControl")
lec <- new("LinStatExpectCovar", ncol(y), ncol(x))

res <- .Call("R_IndependenceTest", y, x, weights,
              lec, varctrl, PACKAGE = "party")
print(res)
wmw <- wilcox.test(xval ~ yf, exact = FALSE, correct = FALSE)
print(wmw)
stopifnot(isequal(res[2], wmw$p.value))

xf <- factor(cut(rnorm(6000), breaks = c(-Inf, -2, 0.5, Inf)))
x <- sapply(levels(xf), function(l) as.numeric(xf == l))
yf <- factor(cut(rnorm(6000), breaks = c(-Inf, -0.5, 1.5, Inf)))
y <- sapply(levels(yf), function(l) as.numeric(yf == l))
weights <- rep(1, nrow(x))

varctrl <- new("VariableControl")
varctrl@teststat <- factor("quad", levels = c("max", "quad"))
print(varctrl)
lec <- new("LinStatExpectCovarMPinv", ncol(x), ncol(y))
lec@rank <- 0
lec@MPinv <- matrix(0, nrow = ncol(x) * ncol(y), ncol = ncol(x) * ncol(y))
lec@svdmem <- new("svd_mem", ncol(x) * ncol(y))

res <- .Call("R_IndependenceTest", x, y, weights,  
              lec, varctrl, PACKAGE = "party")
print(res)
chis <- chisq.test(table(xf, yf), correct = FALSE)
print(chis)
stopifnot(isequal(round(res[2],3), round(chis$p.value,3)))


### Multiple Variables
gtctrl <- new("GlobalTestControl")
tlev <- levels(gtctrl@testtype)   
gtctrl@testtype <- factor("Univariate", levels = tlev)

mydata <- data.frame(y = gl(2, 50), x1 = rnorm(100),
                    x2 = rnorm(100), x3 = rnorm(100))
inp <- initVariableFrame(mydata[,c("x1", "x2", "x3"),drop = FALSE], 
    trafo = function(data) ptrafo(data, numeric_trafo = rank))
resp <- initVariableFrame(mydata[,"y",drop = FALSE], trafo = NULL, response = TRUE)
ls <- new("LearningSample", inputs = inp, responses = resp,
          weights = rep(1, inp@nobs), nobs = nrow(mydata), 
          ninputs = inp@ninputs)
tm <- ctree_memory(ls)
varctrl <- new("VariableControl")
pvals <- .Call("R_GlobalTest", ls, ls@weights, tm, varctrl, gtctrl, PACKAGE = "party")[[2]]
wpvals <- rep(0, 3)
wpvals[1] <- wilcox.test(x1 ~ y, data = mydata,
                         correct = FALSE, exact = FALSE)$p.value
wpvals[2] <- wilcox.test(x2 ~ y, data = mydata, 
                         correct = FALSE, exact = FALSE)$p.value
wpvals[3] <- wilcox.test(x3 ~ y, data = mydata, 
                         correct = FALSE, exact = FALSE)$p.value
stopifnot(isequal(wpvals, 1 - pvals))

varctrl <- new("VariableControl")
gtctrl@testtype <- factor("MonteCarlo", levels = tlev)
gtctrl@nresample <- as.integer(19999)
inp <- initVariableFrame(mydata[,"x1",drop = FALSE], trafo = function(data)
    ptrafo(data, numeric_trafo = rank))
resp <- initVariableFrame(mydata[,"y",drop = FALSE], trafo = NULL, response = TRUE)
ls <- new("LearningSample", inputs = inp, responses = resp,
          weights = rep(1, inp@nobs), nobs = nrow(mydata), 
          ninputs = as.integer(1))
tm <- ctree_memory(ls)
pvals <- .Call("R_GlobalTest", ls, ls@weights, tm, varctrl, gtctrl, PACKAGE = "party")[[2]]
stopifnot(abs((1 - pvals) - wilcox.test(x1 ~ y, data = mydata, 
    exact = TRUE)$p.value) < 1e-2)
