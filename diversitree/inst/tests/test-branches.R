## Test the branches calculations, with both the R and C++ interfaces.
source("helper-diversitree.R")

context("Branches")

pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(4)
phy <- tree.bisse(pars, max.t=30, x0=0)
states <- phy$tip.state
pars.use <- pars + runif(6, 0, 0.01)

lik <- make.bisse(phy, states)
cmp <- attr(lik(pars.use, intermediates=TRUE), "intermediates")

control <- diversitree:::check.control.ode(list())
info    <- diversitree:::make.info.bisse(NULL)
branches <- diversitree:::make.branches.dtlik(info, control)

cache <- diversitree:::make.cache.bisse(phy, states)

initial.conditions <- diversitree:::initial.conditions.bisse
preset <- NULL

all.branches.matrix <- diversitree:::all.branches.matrix

ans <- all.branches.matrix(pars.use, cache, initial.conditions,
                           branches, preset)

v <- intersect(names(ans), names(cmp))
expect_that(ans[v],
            equals(cmp[v]))
