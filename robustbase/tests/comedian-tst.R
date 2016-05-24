require(robustbase)
source(system.file("test-tools-1.R", package="Matrix", mustWork=TRUE))
showProc.time()

data(hbk); hbk.x <- data.matrix(hbk[, 1:3])

covComed(hbk.x)
covComed(hbk.x, n.iter=4)
showProc.time()

data(radarImage)
covComed(radarImage)
covComed(radarImage[,3:5], n.iter = 5)
showProc.time()

data(bushfire) ; covComed(bushfire)
data(heart);     covComed(heart[, 1:2])
data(starsCYG);  covComed(starsCYG)
data(stackloss); covComed(stack.x)
showProc.time()

if(!robustbase:::doExtras()) quit()

## if ( doExtras ) -----------------------------------------------------------------
## ==============



i.rr <- c("raw.cov", "raw.center", "cov", "center")
n <- 1024 ; p <- 7
set.seed(47)
showSys.time(
    rX <- replicate(100, covComed(matrix(rnorm(n*p), n,p))[i.rr],
                    simplify=FALSE))
## Computing simulation-average (cov / center)  <==> looking at Bias
## _FIXME_ Really look at "MSE = Var + Bias^2" -- or something like
## "simulation-average Squared Error or other Loss"
C0 <- Reduce("+", lapply(rX, `[[`, "raw.cov")) / length(rX)
C. <- Reduce("+", lapply(rX, `[[`,     "cov")) / length(rX)
round(1000 * C0)
round(1000 * C.)
assert.EQ(C0, diag(p), tol= 0.04, giveRE=TRUE) #-> 0.02805
assert.EQ(C., diag(p), tol= 0.09, giveRE=TRUE) #-> 0.06475
## Hmm.. raw.cov is better than cov ??
c00 <- Reduce("+", lapply(rX, `[[`, "raw.center")) / length(rX)
c0  <- Reduce("+", lapply(rX, `[[`,     "center")) / length(rX)
stopifnot(print(sqrt(mean( (c00 - rep(0, p))^2 ))) < 0.005)# 0.004188
stopifnot(print(sqrt(mean( (c0  - rep(0, p))^2 ))) < 0.005)# 0.003434

n <- 4096 ; p <- 11
set.seed(17)
showSys.time(
    r4 <- replicate(64, covComed(matrix(10+rnorm(n*p), n,p))[i.rr],
                    simplify=FALSE))
C0 <- Reduce("+", lapply(r4, `[[`, "raw.cov")) / length(r4)
C. <- Reduce("+", lapply(r4, `[[`,     "cov")) / length(r4)
round(1000 * C0)
round(1000 * C.)
assert.EQ(C0, diag(p), tol = 0.025, giveRE=TRUE) # 0.0162
assert.EQ(C., diag(p), tol = 0.06 , giveRE=TRUE) # 0.0486
## Again... raw.cov better than cov ??
c00 <- Reduce("+", lapply(r4, `[[`, "raw.center")) / length(r4)
c0  <- Reduce("+", lapply(r4, `[[`,     "center")) / length(r4)
assert.EQ(c00, rep(10, p), tol = 2e-4, giveRE=TRUE)# 7.97267e-05 = "raw" is better ?
assert.EQ(c0 , rep(10, p), tol = 2e-4, giveRE=TRUE)# 0.0001036

