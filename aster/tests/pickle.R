
library(aster)

options(digits=4) # avoid rounding differences

data(radish)

pred <- c(0,1,2)
fam <- c(1,3,2)

### need object of type aster to supply to penmlogl and pickle

aout <- aster(resp ~ varb + fit : (Site * Region + Block + Pop),
    pred, fam, varb, id, root, data = radish)

### model matrices for fixed and random effects

modmat.fix <- model.matrix(resp ~ varb + fit : (Site * Region),
    data = radish)
modmat.blk <- model.matrix(resp ~ 0 + fit:Block, data = radish)
modmat.pop <- model.matrix(resp ~ 0 + fit:Pop, data = radish)

rownames(modmat.fix) <- NULL
rownames(modmat.blk) <- NULL
rownames(modmat.pop) <- NULL

idrop <- match(aout$dropped, colnames(modmat.fix))
idrop <- idrop[! is.na(idrop)]
modmat.fix <- modmat.fix[ , - idrop]

nfix <- ncol(modmat.fix)
nblk <- ncol(modmat.blk)
npop <- ncol(modmat.pop)

### try penmlogl

sigma.start <- c(1, 1)

alpha.start <- aout$coefficients[match(colnames(modmat.fix),
    names(aout$coefficients))]
parm.start <- c(alpha.start, rep(0, nblk + npop))

tout <- trust(objfun = penmlogl, parm.start, rinit = 1, rmax = 10,
    sigma = sigma.start, fixed = modmat.fix,
    random = list(modmat.blk, modmat.pop), obj = aout)

eff <- tout$argument * tout$scale
eff.blk <- eff[seq(nfix + 1, nfix + nblk)]
eff.pop <- eff[seq(nfix + nblk + 1, nfix + nblk + npop)]
sigma.crude <- sqrt(c(var(eff.blk), var(eff.pop)))

pout <- pickle(sigma.crude, tout$argument, fixed = modmat.fix,
    random = list(modmat.blk, modmat.pop), obj = aout)

# try optim with method = "Nelder-Mead" and pickle

cache <- new.env(parent = emptyenv())
oout <- optim(sigma.crude, pickle, parm = tout$argument,
    fixed = modmat.fix, random = list(modmat.blk, modmat.pop),
    obj = aout, cache = cache)
oout$convergence == 0
oout$par

sigma.mle <- oout$par
tout <- trust(objfun = penmlogl, cache$parm, rinit = 1, rmax = 10,
    sigma = sigma.mle, fixed = modmat.fix,
    random = list(modmat.blk, modmat.pop), obj = aout)
stopifnot(tout$converged)
parm.mle <- tout$argument

# try pickle1

zwz <- makezwz(sigma.mle, parm.mle, fixed = modmat.fix,
    random = list(modmat.blk, modmat.pop), obj = aout)

pout <- pickle1(sigma.mle, parm.mle, fixed = modmat.fix,
    random = list(modmat.blk, modmat.pop), obj = aout, zwz = zwz, deriv = 1)

foo <- function(sigma) {
    pout <- pickle1(sigma.mle, parm.mle, fixed = modmat.fix,
        random = list(modmat.blk, modmat.pop), obj = aout, zwz = zwz,
        deriv = 1)
    result <- pout$value
    attr(result, "gradient") <- pout$gradient
    return(result)
}

nout <- nlm(foo, sigma.mle)
nout$code
nout$iterations

nrand <- c(nblk, npop)
qux <- function(parm, sigma, zwz) {
    foo <- function(alphaceesigma) pickle3(alphaceesigma, fixed = modmat.fix,
        random = list(modmat.blk, modmat.pop), obj = aout, zwz = zwz, deriv = 2)
    sigma <- as.vector(sigma)
    parm <- as.vector(parm)
    zwz <- as.matrix(zwz)
    iter <- NULL
    repeat {
        tout <- trust(foo, c(parm, sigma), rinit = 1, rmax = 10)
        stopifnot(tout$converged)
        iter <- c(iter, tout$iterations)
        sigma.old <- sigma
        sigma <- tout$argument[nfix + sum(nrand) + seq(along = nrand)]
        parm <- tout$argument[seq(1, nfix + sum(nrand))]
        zwz <- makezwz(sigma, parm, fixed = modmat.fix,
            random = list(modmat.blk, modmat.pop), obj = aout)
        # cat("iteration", length(iter), ":",
        #     all.equal(sigma, sigma.old, check.attributes = FALSE), "\n")
        if (isTRUE(all.equal(sigma, sigma.old))) break
    }
    return(list(sigma = sigma, parm = parm, zwz = zwz, iterations = iter))
}

qout <- qux(parm.mle, sigma.mle, zwz)
qout$iterations
qout$sigma

sigma.mle <- qout$sigma
parm.mle <- qout$parm
zwz.mle <- qout$zwz
alpha.mle <- parm.mle[1:nfix]
c.mle <- parm.mle[-(1:nfix)]
a.mle <- rep(sigma.mle, times = nrand)
b.mle <- a.mle * c.mle

# use optim to get hessian of q(alpha, sigma)
# note: nice analytic formula in inst/doc/re.pdf doesn't work
# with inexact computer arithmetic due to catastrophic cancellation

alphasigma.mle <- c(alpha.mle, sigma.mle)

objfun <- function(alphasigma) pickle2(alphasigma, parm = c.mle,
    fixed = modmat.fix, random = list(modmat.blk, modmat.pop),
    obj = aout, zwz = zwz.mle)$value
gradfun <- function(alphasigma) pickle2(alphasigma, parm = c.mle,
    fixed = modmat.fix, random = list(modmat.blk, modmat.pop),
    obj = aout, zwz = zwz.mle, deriv = 1)$gradient

oout <- optim(alphasigma.mle, objfun, gradfun, method = "BFGS", hessian = TRUE)
oout$convergence
hessian.mle <- oout$hessian

