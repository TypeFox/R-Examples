### R code from vignette source '/home/es/Packages/NMOF/reports/NMOFdist/NMOFdist.Rnw'

###################################################
### code chunk number 1: NMOFdist.Rnw:77-78
###################################################
options(continue = " ", digits = 3, width = 65)


###################################################
### code chunk number 2: NMOFdist.Rnw:111-113 (eval = FALSE)
###################################################
## install.packages("NMOF") ## CRAN
## install.packages("NMOF", repos = "http://enricoschumann.net/R")


###################################################
### code chunk number 3: NMOFdist.Rnw:117-120
###################################################
require("NMOF")
set.seed(1122344)
nC <- 4L ## the number of cores to be used


###################################################
### code chunk number 4: NMOFdist.Rnw:126-128
###################################################
require("rbenchmark")
require("parallel")


###################################################
### code chunk number 5: NMOFdist.Rnw:134-136 (eval = FALSE)
###################################################
## whereToLook <- system.file("NMOFex/NMOFdist.R", package = "NMOF")
## file.show(whereToLook, title = "NMOF examples")


###################################################
### code chunk number 6: NMOFdist.Rnw:148-157
###################################################
testFun <- function(ignore, delay) {
    Sys.sleep(delay)
    1
}
 
delay <- 0.05     ## running time of function
n <- 8            ## how many calls per lapply
repl <- 10        ## how many restarts
sq <- seq_len(n)


###################################################
### code chunk number 7: NMOFdist.Rnw:161-168
###################################################
cl <- makeCluster(c(rep("localhost", nC)), type = "SOCK")
benchmark(lapply(sq, testFun, delay),                ## serial
          mclapply(sq, testFun, delay),              ## formerly 'multicore'
          clusterApply(cl, sq, testFun, delay),      ## formerly 'snow'
          columns = c("test", "elapsed", "relative"),
          order = "relative", replications = repl)
stopCluster(cl)


###################################################
### code chunk number 8: NMOFdist.Rnw:176-180
###################################################
OF <- function(b, X, y) {
    temp <- X %*% b - y
    sum(temp^2)
}


###################################################
### code chunk number 9: NMOFdist.Rnw:187-191
###################################################
ncol <- 10; nrow <- 200
X <- array(rnorm(nrow * ncol), dim = c(nrow, ncol))
y <- rnorm(nrow)
b <- rnorm(ncol)


###################################################
### code chunk number 10: NMOFdist.Rnw:195-196
###################################################
OF(b, X, y)


###################################################
### code chunk number 11: NMOFdist.Rnw:201-206
###################################################
n <- 50            ## how many calls per lapply
sq <- seq_len(n)
lP <- vector("list", length = n)
for (i in sq)
    lP[[i]] <- b


###################################################
### code chunk number 12: NMOFdist.Rnw:210-234
###################################################
snow_with_copying <- expression({
    ignore1 <- clusterApply(cl, lP, OF, X, y)  
})

lapply <- expression({
    ignore2 <- lapply(lP, OF, X, y)  
})

snow_without_copying <- expression({
    OF1 <- function(b) {
        temp <- X %*% b - y
        sum(temp^2)
    }
    ignore3 <- clusterApply(cl, lP, OF1)
})

cl <- makeCluster(rep("localhost", nC), type = "SOCK")
clusterExport(cl, list("X", "y"))
benchmark(lapply,            
          snow_with_copying, 
          snow_without_copying,
          columns = c("test", "elapsed", "relative"),
          order = "relative", replications = 50)
stopCluster(cl)


###################################################
### code chunk number 13: NMOFdist.Rnw:241-243
###################################################
all.equal(ignore1, ignore2)
all.equal(ignore2, ignore3)


###################################################
### code chunk number 14: NMOFdist.Rnw:252-278
###################################################
testFun <- function(x) {
    Sys.sleep(0.1)
    cos(1/x^2)
}
with_loop <- expression(
    sol1 <- bracketing(testFun,
                       interval = c(0.3, 0.9),
                       n = 100L))

with_multicore <- expression(
    sol2 <- bracketing(testFun,
                       interval = c(0.3, 0.9),
                       n = 100L,
                       method = "multicore", 
                       mc.control = list(mc.cores = nC)))

with_snow  <- expression(
    sol3 <- bracketing(testFun,
                       interval = c(0.3, 0.9),
                       n = 100L, method = "snow", cl = nC))

benchmark(with_loop, 
          with_multicore, 
          with_snow,
          columns = c("test", "elapsed", "relative"),
          order = "relative", replications = 1)


###################################################
### code chunk number 15: NMOFdist.Rnw:282-284
###################################################
all.equal(sol1, sol2)
all.equal(sol1, sol3)


###################################################
### code chunk number 16: NMOFdist.Rnw:294-297
###################################################
ncol <- 20
nrow <- 1000
P <- array(rnorm(nrow * ncol), dim = c(nrow, ncol))


###################################################
### code chunk number 17: NMOFdist.Rnw:304-306
###################################################
fun <- function (x, h)
    sort(x, partial = h)[h]


###################################################
### code chunk number 18: NMOFdist.Rnw:310-312
###################################################
h <- 5L
fun(P[ ,1L], h)


###################################################
### code chunk number 19: NMOFdist.Rnw:320-327
###################################################
loopfun <- function(x, f, ...) {
    ns <- ncol(x)
    fv <- numeric(ns)
    for (i in seq_len(ns))
        fv[i] <- f(x[ ,i], ...)
    fv
}


###################################################
### code chunk number 20: NMOFdist.Rnw:332-333
###################################################
loopresult <- loopfun(P, fun, h)


###################################################
### code chunk number 21: NMOFdist.Rnw:343-351
###################################################
mat2list <- function(x) {
    nx <- ncol(x)
    listP <- vector(mode = "list", length = nx)
    for (s in seq_len(nx))
        listP[[s]] <- P[ ,s]
    listP
}
listP <- mat2list(P)


###################################################
### code chunk number 22: NMOFdist.Rnw:358-362
###################################################
cl <- makeCluster(c(rep("localhost", nC)), type = "SOCK")
snowresult <- unlist(clusterApply(cl, listP, fun, h))
stopCluster(cl)
all.equal(loopresult, snowresult)


###################################################
### code chunk number 23: NMOFdist.Rnw:367-373
###################################################
cl <- makeCluster(c(rep("localhost", nC)), type = "SOCK")
benchmark(clusterApply(cl, listP, fun, h),
          loopfun(P, fun, h),
          columns = c("test", "elapsed", "relative"),
          order = "relative", replications = 100)
stopCluster(cl)


###################################################
### code chunk number 24: NMOFdist.Rnw:383-388
###################################################
ncol <- 100
nrow <- 1000
P <- array(rnorm(nrow * ncol), dim = c(nrow, ncol))

system.time(for (i in seq_len(10000L)) fun(P[ ,1L], 10L))


###################################################
### code chunk number 25: NMOFdist.Rnw:392-393
###################################################
d <- round(ncol/nC) ## nC is the number of cores


###################################################
### code chunk number 26: NMOFdist.Rnw:396-399
###################################################
listP <- vector(mode = "list", length = nC)
for (s in seq_len(nC))
    listP[[s]] <- P[ ,(d*s-d+1):min(ncol, d*s)]


###################################################
### code chunk number 27: NMOFdist.Rnw:403-410
###################################################
cl <- makeCluster(c(rep("localhost", nC)), type = "SOCK")
benchmark(parallel.result <- clusterApply(cl, listP, loopfun, fun, h),
          loop.result <- loopfun(P, fun, h),
          columns = c("test", "elapsed", "relative"),
          order = "relative", replications = 100)
stopCluster(cl)
all.equal(loop.result, unlist(parallel.result))


###################################################
### code chunk number 28: NMOFdist.Rnw:427-449
###################################################
OF <- function(x, y) {
    Sys.sleep(0.001)
    sum(x != y)
}
size <- 20L            ## the length of the string
y <- runif(size) > 0.5 ## the true solution
with_loop <- list(nB = size, nP = 200L, nG = 50L, prob = 0.002,
                  printBar = FALSE, printDetail = FALSE,
                  methodOF = "loop")
with_snow <- list(nB = size, nP = 200L, nG = 50L, prob = 0.002,
                  printBar = FALSE, printDetail = FALSE,
                  methodOF = "snow", cl = nC)
with_multicore <- list(nB = size, nP = 200L, nG = 50L, prob = 0.002,
                       printBar = FALSE, printDetail = FALSE,
                       methodOF = "multicore")

benchmark(GAopt(OF, algo = with_loop, y = y),
          GAopt (OF, algo = with_snow, y = y),
          GAopt(OF, algo = with_multicore, y = y),
          columns = c("test", "elapsed", "relative"),
          order = "relative", replications = 1)



###################################################
### code chunk number 29: NMOFdist.Rnw:457-463
###################################################
with_multicore$mc.control <- list(mc.cores = 1L)
## system.time(GAopt(OF, algo = with_multicore, y = y))
benchmark(GAopt(OF, algo = with_loop, y = y),
          GAopt(OF, algo = with_multicore, y = y),
          columns = c("test", "elapsed", "relative"),
          order = "relative", replications = 1)


###################################################
### code chunk number 30: NMOFdist.Rnw:467-483
###################################################
OF <- function(x, y) {
    Sys.sleep(0.01)
    sum(x != y)
}
size <- 10L; y <- runif(size) > 0.5
algo <- list(nB = size, nP = 20L, nG = 100L, prob = 0.002,
             printBar = FALSE, methodOF = "loop")
t1 <- system.time(sol <- GAopt(OF, algo = algo, y = y))
all.equal(sol$xbest, y)
all.equal(sol$OFvalue, 0)

algo <- list(nB = size, nP = 20L, nG = 100L, prob = 0.002,
             printBar = FALSE, methodOF = "snow", cl = nC)
t2 <- system.time(sol <- GAopt(OF, algo = algo, y = y))
all.equal(sol$xbest, y)
all.equal(sol$OFvalue, 0)


###################################################
### code chunk number 31: NMOFdist.Rnw:487-488
###################################################
round(t1[[3L]]/t2[[3L]],1)


###################################################
### code chunk number 32: NMOFdist.Rnw:492-510
###################################################
OF <- function(x, y, k) {
    Sys.sleep(0.01)
    sum(x != y) + k
}
size <- 10L; y <- runif(size) > 0.5; k <- 10
algo <- list(nB = size, nP = 20L, nG = 100L, prob = 0.002,
             printBar = FALSE, printDetail = FALSE,
             methodOF = "loop")
t1 <- system.time(sol <- GAopt(OF, algo = algo, y = y, k = k))
all.equal(sol$xbest, y)
all.equal(sol$OFvalue, k)

algo <- list(nB = size, nP = 20L, nG = 100L, prob = 0.002,
             printBar = FALSE, printDetail = FALSE,
             methodOF = "snow", cl = nC)
t2 <- system.time(sol <- GAopt(OF, algo = algo, y = y, k = k))
all.equal(sol$xbest,y)
all.equal(sol$OFvalue, k)


###################################################
### code chunk number 33: NMOFdist.Rnw:516-545
###################################################
testFun  <- function(x) {
    Sys.sleep(0.1)
    x[1L] + x[2L]^2
}
lower <- 1:2; upper <- 5; n <- 10
with_loop <- expression(
    sol1 <- gridSearch(fun = testFun,
                       lower = lower, upper = upper,
                       n = n, printDetail = FALSE))

with_multicore <- expression(
    sol2 <- gridSearch(fun = testFun,
                       lower = lower, upper = upper,
                       n = n, printDetail = FALSE,
                       method = "multicore"))

with_snow <- expression(
    sol3 <- gridSearch(fun = testFun,
                       lower = lower, upper = upper,
                       n = n, printDetail = FALSE,
                       method = "snow",
                       cl = nC))

benchmark(with_loop, with_multicore, with_snow,
          columns = c("test", "elapsed", "relative"),
          order = "relative", replications = 1)
all.equal(sol1, sol2)
all.equal(sol1, sol3)
all.equal(sol3$minlevels, 1:2)


###################################################
### code chunk number 34: NMOFdist.Rnw:550-569
###################################################
testFun  <- function(x, k) {
    Sys.sleep(0.1)
    x[1L] + x[2L]^2 + k
}
lower <- 1:2; upper <- 5; n <- 5; k <- 1
sol1 <- gridSearch(fun = testFun, k = k,
                   lower = lower, upper = upper,
                   n = n, printDetail = FALSE)
sol2 <- gridSearch(fun = testFun,k = k,
                   lower = lower, upper = upper,
                   n = n, printDetail = FALSE,
                   method = "multicore")
sol3 <- gridSearch(fun = testFun,k = k,
                   lower = lower, upper = upper,
                   n = n, printDetail = FALSE,
                   method = "snow", cl = nC)
all.equal(sol1, sol2)
all.equal(sol1, sol3)
all.equal(sol3$minlevels, 1:2)


###################################################
### code chunk number 35: NMOFdist.Rnw:575-594
###################################################
testFun  <- function(x) {
    Sys.sleep(0.1)
    x[1L] + x[2L] + runif(1)
}
lower <- 1:2; upper <- 5; n <- 3
set.seed(5)
sol2 <- gridSearch(fun = testFun,
                   lower = lower, upper = upper,
                   n = n, printDetail = FALSE,
                   method = "multicore",
                   mc.control = list(mc.set.seed = FALSE))
temp <- sol2$values
set.seed(5)
sol2 <- gridSearch(fun = testFun,
                   lower = lower, upper = upper,
                   n = n, printDetail = FALSE,
                   method = "multicore",
                   mc.control = list(mc.set.seed = FALSE))
all.equal(sol2$values, temp)


###################################################
### code chunk number 36: NMOFdist.Rnw:599-615
###################################################
cl <- makeCluster(c(rep("localhost", nC)), type = "SOCK")
clusterSetRNGStream(cl, 2222)
sol3 <- gridSearch(fun = testFun, lower = lower, upper = upper,
                   n = n, printDetail = FALSE,
                   method = "snow", cl = cl)
stopCluster(cl)
temp <- sol3$values

## ... and again
cl <- makeCluster(c(rep("localhost", nC)), type = "SOCK")
clusterSetRNGStream (cl, 2222)
sol3 <- gridSearch(fun = testFun, lower = lower, upper = upper,
                   n = n, printDetail = FALSE,
                   method = "snow", cl = cl)
stopCluster(cl)
all.equal(sol3$values, temp)


###################################################
### code chunk number 37: NMOFdist.Rnw:623-635
###################################################
xTRUE <- runif(5L)
data <- list(xTRUE = xTRUE,  ## the TRUE solution
             step = 0.02     ## step size for neighbourhood
             )
OF <- function(x, data)
    max(abs(x - data$xTRUE))
neighbour <- function(x, data)
    x + runif(length(data$xTRUE))*data$step - data$step/2
x0 <- runif(5L)              ## a random starting solution
algo <- list(q = 0.05, nS = 200L, nT = 10L,
             neighbour = neighbour, x0 = x0,
             printBar = FALSE, printDetail = FALSE)


###################################################
### code chunk number 38: NMOFdist.Rnw:638-658
###################################################
with_loop <- expression(
    sols1 <- restartOpt(fun = TAopt, n = 100L,
                        OF = OF, algo = algo, data = data))

with_multicore <- expression(
    sols2 <- restartOpt(fun = TAopt, n = 100L,
                        OF = OF, algo = algo, data = data,
                        method = "multicore"))

with_snow <- expression(
    sols3 <- restartOpt(fun = TAopt, n = 100L,
                        OF = OF, algo = algo, data = data,
                        method = "snow", cl = nC))

benchmark(with_loop, with_multicore, with_snow,
          columns = c("test", "elapsed", "relative"),
          order = "relative", replications = 1)
all.equal(length(sols1), 100L)
all.equal(length(sols2), 100L)
all.equal(length(sols3), 100L)


###################################################
### code chunk number 39: NMOFdist.Rnw:685-686
###################################################
toLatex(sessionInfo())


