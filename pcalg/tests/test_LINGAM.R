library(pcalg)
## source("/u/kalischm/research/packages/LINGAM/R/lingamFuns.R")

##--> showProc.time(), assertError(), relErrV(), ...
R.home(); sessionInfo() # helping package maintainers to debug ...
.libPaths()
packageDescription("pcalg")
packageDescription("Matrix")
cat("doExtras:", (doExtras <- pcalg:::doExtras()), "\n")


##################################################
## Exp 1
##################################################
set.seed(123)
n <- 500
eps1 <- sign(rnorm(n)) * sqrt(abs(rnorm(n)))
eps2 <- runif(n) - 0.5 #  ~ U[-1/2, 1/2]

X <- cbind(A = eps1 + 0.9*eps2,
           B =            eps2)

## x1 <- x2
## adjacency matrix:
## 0 0
## 1 0
(trueDAG <- rbind(c(0,0),
                  c(1,0)))

estDAG <- LINGAM(X, verbose = TRUE)

stopifnot(as.integer(estDAG$ Adj) == trueDAG,
          all.equal (estDAG$ B, cbind(0, c(0.878188262685122, 0))))

if(doExtras) {
    ## using pcalg
    n <- nrow(X)
    V <- LETTERS[1:ncol(X)] # labels aka node names

    ## estimate CPDAG
    pc.fit <- pc(suffStat = list(C = cor(X), n = n),
                 indepTest = gaussCItest, ## indep.test: partial correlations
                 alpha = 0.01, labels = colnames(X))
    if (require(Rgraphviz)) {
        plot(pc.fit, main = "Estimated CPDAG")
    }
}

##################################################
## Exp 2
##################################################
set.seed(123)
n <- 500
eps1 <- sign(rnorm(n)) * sqrt(abs(rnorm(n)))
eps2 <- runif(n) - 0.5
eps3 <- sign(rnorm(n)) * abs(rnorm(n))^(1/3)
eps4 <- rnorm(n)^2

x1 <- eps1 + 0.9*eps2
x2 <-            eps2
x3 <-        0.8*eps2 + eps3
x4 <- -0.9*x3 - x1 + eps4

X <- cbind(U = x1, V = x2, W = x3, Y = x4)

trueDAG <- cbind(c(0,1,0,0),c(0,0,0,0),c(0,1,0,0),c(1,0,1,0))
## x4 <- x3 <- x2 -> x1 -> x4
## adjacency matrix:
## 0 0 0 1
## 1 0 1 0
## 0 0 0 1
## 0 0 0 0

estDAG <- LINGAM(X, verbose = TRUE)

B.est <- rbind(c(0, 0.986119553, 0, 0),
               c(0, 0, 0, 0),
               c(0, 0.89198226, 0, 0),
               c(-0.987301824, 0, -0.890961952, 0))

stopifnot(as.integer(estDAG$Adj) == trueDAG,
          all.equal(estDAG$B, B.est, tol=1e-9))

if(doExtras) {
    ## using pcalg
    n <- nrow(X)
    V <- colnames(X) # labels aka node names

    ## estimate CPDAG
    pc.fit <- pc(suffStat = list(C = cor(X), n = n),
                 indepTest = gaussCItest, ## indep.test: partial correlations
                 alpha=0.01, labels = V, verbose = FALSE)
    if (require(Rgraphviz)) {
        plot(pc.fit, main = "Estimated CPDAG")
    }
}

## if(!doExtras && !interactive()) quit("no")

### More tests for higher dimensions

### p = 8 -------- Example 3 -----------

set.seed(127)
n <- 2000
x1 <- eps1 <- sign(rnorm(n)) * sqrt(abs(rnorm(n)))
x2 <- eps2 <- runif(n) - 0.5
x3 <- eps3 <- sign(rnorm(n)) * abs(rnorm(n))^(1/3)
x4 <- eps4 <- rnorm(n)^2
Z <- rnorm(n); eps5 <- sign(Z) * sqrt(abs(Z))
Z <- rnorm(n); eps6 <- sign(Z) * sqrt(abs(Z))
Z <- rnorm(n); eps7 <- sign(Z) * sqrt(abs(Z))
Z <- rnorm(n); eps8 <- sign(Z) * sqrt(abs(Z))

x5 <- 7/8*x1 - 3/4*x2 + 3/4*x3 + eps5
x6 <- 0.8*x4                   + eps6
x7 <- 3/4*x5 - 7/8*x6          + eps7
x8 <- .9*x7                    + eps8

X <- cbind(x1,x2,x3,x4,x5,x6,x7,x8, deparse.level = 2)
## (x1, x2, x3) -> x5 -> x7 <- x6 <- x4;  x7 -> x8
## adjacency matrix:
##     1 2 3 4 5 6 7 8
## x1  . . . . 1 . . .
## x2  . . . . 1 . . .
## x3  . . . . 1 . . .
## x4  . . . . . 1 . .
## x5  . . . . . . 1 .
## x6  . . . . . . 1 .
## x7  . . . . . . . 1
## x8  . . . . . . . .

## true DAG :
. <- 0
trDAG3 <- rbind(
    c(., ., ., ., 1, ., ., .),
    c(., ., ., ., 1, ., ., .),
    c(., ., ., ., 1, ., ., .),
    c(., ., ., ., ., 1, ., .),
    c(., ., ., ., ., ., 1, .),
    c(., ., ., ., ., ., 1, .),
    c(., ., ., ., ., ., ., 1),
    c(., ., ., ., ., ., ., .))

estB.3 <- rbind(
    c(., ., ., ., ., ., ., .),
    c(., ., ., ., ., ., ., .),
    c(., ., ., ., ., ., ., .),
    c(., ., ., ., ., ., ., .),
    c(.831899433, -.737954104, 0.725137273, ., ., ., ., .),
    c(., ., ., 0.788185348, ., ., ., .),
    c(., ., ., ., 0.774490692, -0.886143314, ., .),
    c(., ., ., ., ., ., 0.900617843, .))

eDAG3 <- LINGAM(X, verbose = TRUE)

stopifnot(trDAG3 == eDAG3$Adj,
          with(eDAG3, all(t(B != 0) == Adj)),
          all.equal(eDAG3$B, estB.3, tol=1e-9))


### p = 10 -------- Example 4 -----------

### using same x1..,x4, and, eps5 .. eps8  as in Ex. 3

Z <- rnorm(n); eps9  <- sign(Z) * abs(Z)^(1/3)
Z <- rnorm(n); eps10 <- sign(Z) * abs(Z)^0.25

x5 <-  7/8*x1 - 3/4*x2          + eps5
x6 <-  0.8*x2 - 7/8*x3          + eps6
x7 <- -7/8*x4                   + eps7
x8 <-  0.9*x2 - 0.8*x5          + eps8
x9 <- -3/4*x6 + 7/8*x7          + eps9
x10 <- 3/4*x6 + 0.5*x8 +0.9*x9  + eps10

X <- cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10, deparse.level = 2)


## true DAG :
. <- 0
trDAG4 <- rbind(
    # 1  2  3  4  5  6  7  8  9  0
    c(., ., ., ., 1, ., ., ., ., .), #  1
    c(., ., ., ., 1, 1, ., 1, ., .), #  2
    c(., ., ., ., ., 1, ., ., ., .), #  3
    c(., ., ., ., ., ., 1, ., ., .), #  4
    c(., ., ., ., ., ., ., 1, ., .), #  5
    c(., ., ., ., ., ., ., ., 1, 1), #  6
    c(., ., ., ., ., ., ., ., 1, .), #  7
    c(., ., ., ., ., ., ., ., ., 1), #  8
    c(., ., ., ., ., ., ., ., ., 1), #  9
    c(., ., ., ., ., ., ., ., ., .)) # 10

estB.4 <- rbind(
    c(., ., ., ., .,  ., ., ., ., .),
    c(., ., ., ., .,  ., ., ., ., .),
    c(., ., ., ., .,  ., ., ., ., .),
    c(., ., ., ., .,  ., ., ., ., .),
    c(0.831899433, -0.737954104, ., ., .,  ., ., ., ., .),
    c(.,            0.687919607, -0.863361084, ., .,  ., ., ., ., .),
    c(.,             .,            .,        -0.878305407, .,  ., ., ., ., .),
    c(.,            0.864092647,   .,          .,        -0.77902333,  ., ., ., ., .),
    c(.,             .,            .,          .,          .,        -0.780116888, 0.929828083, ., ., .),
    c(.,             .,            .,          .,          .,         0.72436897,   .,     0.502210828, 0.913644804, .))


eDAG4 <- LINGAM(X, verbose = TRUE)

stopifnot(trDAG4 == eDAG4$Adj,
          with(eDAG4, all(t(B != 0) == Adj)),
          all.equal(eDAG4$B, estB.4, tol=1e-9))
