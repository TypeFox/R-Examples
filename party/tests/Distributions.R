
set.seed(290875)
library("party")
if (!require("mvtnorm"))
    stop("cannot load package mvtnorm")


### get rid of the NAMESPACE
attach(asNamespace("party"))

### 
###
###    Regression tests for conditional distributions
###    
###    functions defined in file `./src/Distributions.c'

### chisq-distribution of quadratic forms
t <- 2.1
df <- 2
storage.mode(t) <- "double"
storage.mode(df) <- "double"
stopifnot(isequal(1 - pchisq(t, df = df), ### P-values!!!
          .Call("R_quadformConditionalPvalue", t, df, PACKAGE = "party")))

stopifnot(isequal(2*pnorm(-t), 
          .Call("R_maxabsConditionalPvalue", t, matrix(1), as.integer(1), 0.0, 0.0, 0.0, PACKAGE = "party")))


maxpts <- 25000
storage.mode(maxpts) <- "integer"
abseps <- 0.0001
releps <- 0
tol <- 1e-10

a <- 1.96
b <- diag(2)

p1 <- .Call("R_maxabsConditionalPvalue", a, b, maxpts, abseps, releps, tol, PACKAGE = "party")
p2 <- pmvnorm(lower = rep(-a,2), upper = rep(a,2), corr = b)
stopifnot(isequal(round(p1, 3), round(1 - p2, 3)))

b <- diag(4)
p1 <- .Call("R_maxabsConditionalPvalue", a, b, maxpts, abseps, releps, tol, PACKAGE = "party")
p2 <- pmvnorm(lower = rep(-a,4), upper = rep(a,4), corr = b)
stopifnot(isequal(round(p1, 3), round(1 - p2, 3)))

b <- diag(4)
b[upper.tri(b)] <- c(0.1, 0.2, 0.3)
b[lower.tri(b)] <- t(b)[lower.tri(b)]
p1 <- .Call("R_maxabsConditionalPvalue", a, b, maxpts, abseps, releps, tol, PACKAGE = "party")
p2 <- pmvnorm(lower = rep(-a,4), upper = rep(a,4), corr = b)
stopifnot(isequal(round(p1, 3), round(1 - p2, 3)))

### Monte-Carlo approximation of P-Values, univariate
mydata = data.frame(y = gl(2, 50), x1 = rnorm(100),  
                    x2 = rnorm(100), x3 = rnorm(100))
inp <- initVariableFrame(mydata[,"x1",drop = FALSE], trafo = function(data)
ptrafo(data, numeric_trafo = rank))
resp <- initVariableFrame(mydata[,"y",drop = FALSE], trafo = NULL, response = TRUE)
ls <- new("LearningSample", inputs = inp, responses = resp,
          weights = rep(1, inp@nobs), nobs = nrow(mydata),
          ninputs = inp@ninputs)
tm <- ctree_memory(ls)
varctrl <- new("VariableControl")
varctrl@teststat <- factor("max", levels = c("max", "quad"))
varctrl@pvalue <- FALSE
gtctrl <- new("GlobalTestControl")
gtctrl@testtype <- factor("MonteCarlo", levels = levels(gtctrl@testtype))
gtctrl@nresample <- as.integer(19999)

pvals <- .Call("R_GlobalTest", ls, ls@weights, tm, varctrl, gtctrl, PACKAGE = "party")
wstat <- abs(qnorm(wilcox.test(x1 ~ y, data = mydata, 
             exact = FALSE, correct = FALSE)$p.value/2))
wpval <- wilcox.test(x1 ~ y, data = mydata, exact = TRUE)$p.value
stopifnot(isequal(wstat, pvals[[1]]))
stopifnot(abs(wpval - (1 - pvals[[2]])) < 0.01)

### Monte-Carlo approximations of P-Values, multiple inputs
mydata = data.frame(y = gl(2, 50), x1 = rnorm(100),  
                    x2 = rnorm(100), x3 = rnorm(100))
inp <- initVariableFrame(mydata[,c("x1", "x2", "x3"),
                                drop = FALSE], trafo = function(data)
ptrafo(data, numeric_trafo = rank))
resp <- initVariableFrame(mydata[,"y",drop = FALSE], trafo = NULL, response = TRUE)
ls <- new("LearningSample", inputs = inp, responses = resp,
          weights = rep(1, inp@nobs), nobs = nrow(mydata),
          ninputs = inp@ninputs)
tm <- ctree_memory(ls)
varctrl <- new("VariableControl")
varctrl@teststat <- factor("max", levels = c("max", "quad"))
varctrl@pvalue <- TRUE
gtctrl <- new("GlobalTestControl")
gtctrl@testtype <- factor("Univariate", levels = levels(gtctrl@testtype))
gtctrl@nresample <- as.integer(19999)

pvals <- .Call("R_GlobalTest", ls, ls@weights, tm, varctrl, gtctrl, PACKAGE = "party")
wstat <- c(abs(qnorm(wilcox.test(x1 ~ y, data = mydata, 
               exact = FALSE, correct = FALSE)$p.value/2)),
           abs(qnorm(wilcox.test(x2 ~ y, data = mydata, 
               exact = FALSE, correct = FALSE)$p.value/2)),
           abs(qnorm(wilcox.test(x3 ~ y, data = mydata, 
               exact = FALSE, correct = FALSE)$p.value/2)))
wpval <- c(wilcox.test(x1 ~ y, data = mydata, 
               exact = FALSE, correct = FALSE)$p.value,
           wilcox.test(x2 ~ y, data = mydata, 
               exact = FALSE, correct = FALSE)$p.value,
           wilcox.test(x3 ~ y, data = mydata, 
               exact = FALSE, correct = FALSE)$p.value)
stopifnot(isequal(wstat, pvals[[1]]))
stopifnot(isequal(wpval, 1 - pvals[[2]]))

### Monte-Carlo approximations of P-Values, min-P approach
gtctrl@testtype <- factor("MonteCarlo", levels = levels(gtctrl@testtype))
gtctrl@nresample <- as.integer(19999)
pvals <- .Call("R_GlobalTest", ls, ls@weights, tm, varctrl, gtctrl, PACKAGE = "party")
stopifnot(isequal(wstat, pvals[[1]]))
stopifnot(all(wpval < (1 - pvals[[2]])))
