library(mcmc)

.morph.unmorph <- mcmc:::.morph.unmorph

###########################################################################
# basic functionality check, can morph.metro run?  Can we change the
# transformation?
set.seed(42)
obj <- morph.metrop(function(x) dt(x, df=3, log=TRUE),
                    100, 100, morph=morph(b=3))
obj <- morph.metrop(obj, morph=morph(b=1))

obj <- morph.metrop(function(x) prod(dt(x, df=3, log=TRUE)),
                    rep(100, 3), 100, morph=morph(p=3, b=1))
obj <- morph.metrop(obj, morph=morph(r=1, p=3, b=1))

all.equal(class(obj), c("mcmc", "morph.metropolis"))

###########################################################################
# check .morph.unmorph
obj <- list(final=10)
outfun <- function(x) x
m <- morph(p=3)
obj <- .morph.unmorph(obj, m, outfun)
all.equal(class(obj), c("mcmc", "morph.metropolis"))
all.equal(sort(names(obj)),
          sort(c("final", "morph", "morph.final", "outfun")))
all.equal(c(obj$final, obj$morph.final), c(m$inverse(10), 10))
all.equal(obj$outfun, outfun)
all.equal(obj$morph, m)
