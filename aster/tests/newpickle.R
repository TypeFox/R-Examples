
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

### try newpickle

sigma.start <- c(1, 1)
alpha.start <- aout$coefficients[match(colnames(modmat.fix),
    names(aout$coefficients))]
cee.start <- rep(0, nblk + npop)
alphaceesigma.start <- c(alpha.start, cee.start, sigma.start)

foo <- newpickle(alphaceesigma.start, fixed = modmat.fix,
    random = list(modmat.blk, modmat.pop), obj = aout)
names(foo)

