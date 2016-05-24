
##----setup2
require(Matrix, quietly=TRUE)
NN <- 6
kk <- 2
pp <- 2
nv1 <- NN*kk+pp
nels1 <- nv1^2
nnz1 <- NN*kk^2 + pp^2 + 2*NN*pp*kk
nnz1LT <- NN*kk*(kk+1)/2 + pp*(pp+1)/2 + pp*NN*kk
Q <- 1000
nv2 <- Q*kk+pp
nels2 <- nv2^2
nnz2 <- Q*kk^2 + pp^2 + 2*Q*kk*pp
nnz2LT <- Q*kk*(kk+1)/2 + pp*(pp+1)/2 + Q*kk*pp



##----pattern1
MM <- as(kronecker(diag(NN),matrix(1,kk,kk)),"lMatrix")
MM <- rBind(MM, Matrix(TRUE,pp,NN*kk))
MM <- cBind(MM, Matrix(TRUE, kk*NN+pp, pp))
print(as(MM,"lgCMatrix"))



##----pattern2
MM <- as(kronecker(matrix(1,kk,kk), diag(NN)),"lMatrix")
MM <- rBind(MM, Matrix(TRUE,pp,NN*kk))
MM <- cBind(MM, Matrix(TRUE, kk*NN+pp, pp))
print(as(MM,"lgCMatrix"))


##----data
data(binary_small)
binary <- binary_small #rename for brevity
str(binary)
N <- length(binary[["Y"]])
k <- NROW(binary[["X"]])
q <- k
nvars <- as.integer(N*k + q)
priors <- list(inv.Sigma = diag(k), ##rWishart(1,k+5,diag(k))[,,1],
               inv.Omega = diag(k))
data.frame(parameter = c("N","k","q"),
           value =  c(N, k, q))



##----startingVals
start <- rnorm(nvars) ## random starting values
f <- binary.f(start, data=binary, priors=priors)
f
df <- binary.grad(start, data=binary, priors=priors)
str(df)
d2f <- binary.hess(start, data=binary, priors=priors)
print(d2f[1:6,1:6], digits=3)


# ---- hessStruct
require(sparseHessianFD, quietly=TRUE)
hs <- Matrix(0, nvars, nvars)
for (i in 1:(N + 1)) {
    ## range of row / col indices of block diagonal
    rng <- ((i-1)*k+1):(k*i)
    hs[rng, rng] <- tril(Matrix(1,k,k)) ## lower triangle
}
hs[N*k + 1:q, 1:(N*k)] <- 1 ## bottom margin
hsNZ <- Matrix.to.Coord(hs)
str(hsNZ)


##----hessStruct2
## require(sparseHessianFD)
## hs <- drop0(tril(binary.hess(start, data=binary, priors=priors)))
## hsNZ <- Matrix.to.Coord(hs)
## str(hsNZ)

##----sparseHessianFD
FD <- sparseHessianFD(start, binary.f, binary.grad,
                          rows=hsNZ[["rows"]], cols=hsNZ[["cols"]],
                          data=binary, priors=priors)



##----usingFD
f <- FD$fn(start)
df <- FD$gr(start)
hess <- FD$hessian(start)



##----hessUpperLeft
print(hess[1:6,1:6], digits=3)
all.equal(hess, d2f, tolerance = 1e-7)


##----trustOptim
require(trustOptim, quietly=TRUE)
opt <- trust.optim(start, fn=FD$fn, gr = FD$gr, hs = FD$hessian,
                   method = "Sparse",
                   control = list(
                       start.trust.radius=5, stop.trust.radius = 1e-7,
                       prec=1e-7, report.precision=1,
                       maxit=500, preconditioner=1,
                       function.scale.factor=-1
                       )
                   )

theta.star <- opt[["solution"]]
hess <- opt[["hessian"]]


##----defPropFuncs
require(sparseMVN, quietly=TRUE)
rmvn.sparse.wrap <- function(n.draws, params) {
    rmvn.sparse(n.draws, params[["mean"]], params[["CH"]], prec=TRUE)
}
dmvn.sparse.wrap <- function(d, params) {
    dmvn.sparse(d, params[["mean"]], params[["CH"]], prec=TRUE)
}



##----propParams
scale <- .96
chol.hess <- Cholesky(-scale*hess)
prop.params <- list(mean = theta.star, CH = chol.hess)

##----parallelSetup
library(doParallel, quietly=TRUE)
run.par <- TRUE
if(run.par) registerDoParallel(cores=2) else registerDoParallel(cores=1)
seed.id <- 123
set.seed(seed.id)



##----proposals
M <- 10000  ## proposal draws
log.c1 <- FD$fn(theta.star)
log.c2 <- dmvn.sparse.wrap(theta.star, prop.params)
draws.m <- as(rmvn.sparse.wrap(M,prop.params),"matrix")
log.post.m <- plyr::aaply(draws.m, 1, FD$fn, .parallel=run.par)
log.prop.m <- dmvn.sparse.wrap(draws.m, params=prop.params)
log.phi <- log.post.m - log.prop.m + log.c2 - log.c1
valid.scale <- all(log.phi <= 0)
stopifnot(valid.scale)


##----sampleGDS_serial
n.draws <- 5  ## total number of draws needed
max.tries <- 100000  ## to keep sample.GDS from running forever
if (!run.par) {
    draws <- sample.GDS(n.draws = n.draws,
                        log.phi = log.phi,
                        post.mode = theta.star,
                        fn.dens.post = FD$fn,
                        fn.dens.prop = dmvn.sparse.wrap,
                        fn.draw.prop = rmvn.sparse.wrap,
                        prop.params = prop.params,
                        report.freq = 1, announce=TRUE)
}

##----sampleGDS_parallel
if (run.par) {
    n.batch <- 10
    batch.size <- ceiling(n.draws / n.batch)
    draws.list <- foreach(i=1:n.batch, .inorder=FALSE) %dopar% sample.GDS(
        n.draws = n.draws,
        log.phi = log.phi,
        post.mode = theta.star,
        fn.dens.post = FD$fn,
        fn.dens.prop = dmvn.sparse.wrap,
        fn.draw.prop = rmvn.sparse.wrap,
        prop.params = prop.params,
        report.freq = 1,
        thread.id = i,
        announce=TRUE,
        seed=as.integer(seed.id*i))
    ## combine results from each batch
    draws <- Reduce(function(x,y) Map(rbind,x,y), draws.list)
}


##----strDraws
str(draws)


##----summary
quants <-  plyr::aaply(draws[["draws"]][,(N*k+1):nvars], 2,
                       quantile, probs=c(.025, .5, .975),
                       .parallel = run.par)
quants




##----LML
if (any(is.na(draws[["counts"]]))) {
    LML <- NA
} else {
    LML <- get.LML(counts=draws$counts,
                   log.phi=log.phi,
                   post.mode=theta.star,
                   fn.dens.post= FD$fn,
                   fn.dens.prop=dmvn.sparse.wrap,
                   prop.params=prop.params)
}
draws[["LML"]] <- LML
draws[["acc.rate"]] <- 1/mean(draws$counts)
cat("Acceptance rate: ",draws[["acc.rate"]],
    "\nLog marginal likelihood: ",draws[["LML"]],"\n")

