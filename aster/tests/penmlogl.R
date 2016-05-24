
library(aster)

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

theta.start <- c(1, 1)

alpha.start <- aout$coefficients[match(colnames(modmat.fix),
    names(aout$coefficients))]
parm.start <- c(alpha.start, rep(0, nblk + npop))

tout <- trust(objfun = penmlogl, parm.start, rinit = 1, rmax = 10,
    sigma = theta.start, fixed = modmat.fix,
    random = list(modmat.blk, modmat.pop), obj = aout)
tout$converged

### crude estimate of variance components

eff.blk <- tout$argument[seq(nfix + 1, nfix + nblk)]
eff.pop <- tout$argument[seq(nfix + nblk + 1, nfix + nblk + npop)]
theta.crude <- sqrt(c(var(eff.blk), var(eff.pop)))

### now for the testing

pout <- penmlogl(tout$argument, theta.crude, fixed = modmat.fix,
    random = list(modmat.blk, modmat.pop), obj = aout)
names(pout)

eff.fix <- tout$argument[seq(1, nfix)]
eff.blk <- tout$argument[seq(nfix + 1, nfix + nblk)]
eff.pop <- tout$argument[seq(nfix + nblk + 1, nfix + nblk + npop)]
eff.blk <- theta.crude[1] * eff.blk
eff.pop <- theta.crude[2] * eff.pop

beta <- c(eff.fix, eff.blk, eff.pop)
scalevec <- rep(c(1, theta.crude), times = c(nfix, nblk, npop))
penaltyvec <- rep(c(0, 1, 1), times = c(nfix, nblk, npop))

modmat <- cbind(modmat.fix, modmat.blk, modmat.pop)
mout <- mlogl(beta, pred, fam, aout$x, aout$root, modmat, deriv = 2)

all.equal(mout$hessian * outer(scalevec, scalevec), pout$mlogl.hessian)
all.equal(tout$argument, pout$argument)
all.equal(scalevec, pout$scale)

all.equal(pout$value,
    mout$value + sum(penaltyvec * pout$argument^2) / 2)
all.equal(pout$gradient,
    mout$gradient * scalevec + penaltyvec * pout$argument)
all.equal(pout$hessian, pout$mlogl.hessian + diag(penaltyvec))

epsilon <- 1e-7
epsilon
mygradient <- 0 * pout$gradient
myhessian <- 0 * pout$hessian
for (i in seq(along = mygradient)) {
    arg <- tout$argument
    arg[i] <- arg[i] + epsilon
    pout.epsilon <- penmlogl(arg, theta.crude, fixed = modmat.fix,
        random = list(modmat.blk, modmat.pop), obj = aout)
    mygradient[i] <- (pout.epsilon$value - pout$value) / epsilon
    myhessian[i, ] <- (pout.epsilon$gradient - pout$gradient) / epsilon
}
all.equal(mygradient, pout$gradient, tol = 100 * epsilon)
all.equal(myhessian, pout$hessian, tol = 100 * epsilon)

### penmlogl2

foo <- pout$argument
alpha2 <- foo[seq(1, nfix)]
parm2 <- foo[- seq(1, nfix)]
pout2 <- penmlogl2(parm2, alpha2, theta.crude, fixed = modmat.fix,
    random = list(modmat.blk, modmat.pop), obj = aout)
identical(pout$value, pout2$value)
identical(pout$argument, pout2$argument)
identical(pout$scale[- seq(1, nfix)], pout2$scale)
identical(pout$mlogl.hessian, pout2$mlogl.hessian)
identical(pout$mlogl.gradient, pout2$mlogl.gradient)
idx <- seq(nfix + 1, nfix + nblk + npop)
identical(pout2$gradient, pout$gradient[idx])
foom <- pout$hessian
foom <- foom[idx, ]
foom <- foom[ , idx]
identical(pout2$hessian, foom)

