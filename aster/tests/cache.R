
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

### crude estimate of variance components

eff.blk <- tout$argument[seq(nfix + 1, nfix + nblk)]
eff.pop <- tout$argument[seq(nfix + nblk + 1, nfix + nblk + npop)]
theta.crude <- sqrt(c(var(eff.blk), var(eff.pop)))

### try optim and pickle

lower <- rep(0, length(theta.crude))
time1 <- system.time(
oout <- optim(theta.crude, pickle,
    parm = tout$argument, fixed = modmat.fix,
    random = list(modmat.blk, modmat.pop), obj = aout)
)

time2 <- system.time(
oout.cache <- optim(theta.crude, pickle,
    parm = tout$argument, fixed = modmat.fix,
    random = list(modmat.blk, modmat.pop), obj = aout,
    cache = new.env())
)

all.equal(oout$par, oout.cache$par, tol = 1e-5)
time2[1] < time1[1]

