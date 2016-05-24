source("helper-diversitree.R")

## Will be useful to have a function for testing tolerance to within
## 1e-7, as that works out to be how accurate most things actually
## are.
equals7 <- function(...)
  equals(..., tolerance=1e-7)

context("ClaSSE")

read.ttn <- function(treefile)
{
    treestr <- readLines(treefile, n=1)
    tree <- read.tree(text=treestr)
    ntips <- Ntip(tree)

    all.states <- read.table(treefile, skip=1)
    tip.states <- all.states[1:ntips, 2]
    names(tip.states) <- all.states[1:ntips, 1]

    return(list(tree = tree, states = tip.states))
}

# sim params = (1.4, 0.6, 0.2, 0.2, 0.5, 0), root state = 0
ttn1 <- read.ttn("bisse-tree.ttn")

# sim params = (1.5, 0.5, 1.0, 0.7, 0.7, 1.5, 1.5), root state = 0
ttn2 <- read.ttn("geosse-tree.ttn")

# sim params (lam = 6.2, p = 0.2, q = 0; symmetric):
# lam111 lam112 lam122 lam211 lam212 lam222    mu1    mu2    q12    q21 
#   4.96   1.24      0      0   1.24   4.96      3      0      0      0
ttn3 <- read.ttn("classe-tree.ttn")

## Different control parameters
control.d <- list(backend="deSolve")
control.g <- list(backend="gslode", compiled=FALSE)
control.G <- list(backend="gslode", compiled=TRUE)
control.i <- list(backend="invalid_backend")

lnL1.b <- make.bisse(ttn1$tree, ttn1$states)
lnL1.0 <- make.classe(ttn1$tree, ttn1$states+1, 2)
lnL1.d <- make.classe(ttn1$tree, ttn1$states+1, 2, control=control.d)
lnL1.g <- make.classe(ttn1$tree, ttn1$states+1, 2, control=control.g)
lnL1.G <- make.classe(ttn1$tree, ttn1$states+1, 2, control=control.G)

lnL1.02 <- constrain(lnL1.0,
                     lambda112~0, lambda122~0, lambda211~0, lambda212~0)
lnL1.d2 <- constrain(lnL1.d,
                     lambda112~0, lambda122~0, lambda211~0, lambda212~0)
lnL1.g2 <- constrain(lnL1.g,
                     lambda112~0, lambda122~0, lambda211~0, lambda212~0)
lnL1.G2 <- constrain(lnL1.G,
                     lambda112~0, lambda122~0, lambda211~0, lambda212~0)

lnL2.geosse <- make.geosse(ttn2$tree, ttn2$states)
lnL2.0 <- make.classe(ttn2$tree, ttn2$states+1, 3)
lnL2.d <- make.classe(ttn2$tree, ttn2$states+1, 3, control=control.d)
lnL2.g <- make.classe(ttn2$tree, ttn2$states+1, 3, control=control.g)
lnL2.G <- make.classe(ttn2$tree, ttn2$states+1, 3, control=control.G)

lnL2.02 <- constrain(lnL2.0, lambda111~0, lambda122~0,
                     lambda133~0, lambda211~0, lambda212~0, lambda213~0,
                     lambda223~0, lambda233~0, lambda311~0, lambda312~0,
                     lambda313~0, lambda322~0, lambda323~0,
                     lambda112~lambda222, lambda113~lambda333, 
                     mu1~0, q23~0, q32~0, q13~mu2, q12~mu3)
lnL2.d2 <- constrain(lnL2.d, lambda111~0, lambda122~0,
                     lambda133~0, lambda211~0, lambda212~0, lambda213~0,
                     lambda223~0, lambda233~0, lambda311~0, lambda312~0,
                     lambda313~0, lambda322~0, lambda323~0,
                     lambda112~lambda222, lambda113~lambda333, 
                     mu1~0, q23~0, q32~0, q13~mu2, q12~mu3)
lnL2.g2 <- constrain(lnL2.g, lambda111~0, lambda122~0,
                     lambda133~0, lambda211~0, lambda212~0, lambda213~0,
                     lambda223~0, lambda233~0, lambda311~0, lambda312~0,
                     lambda313~0, lambda322~0, lambda323~0,
                     lambda112~lambda222, lambda113~lambda333, 
                     mu1~0, q23~0, q32~0, q13~mu2, q12~mu3)
lnL2.G2 <- constrain(lnL2.G, lambda111~0, lambda122~0,
                     lambda133~0, lambda211~0, lambda212~0, lambda213~0,
                     lambda223~0, lambda233~0, lambda311~0, lambda312~0,
                     lambda313~0, lambda322~0, lambda323~0,
                     lambda112~lambda222, lambda113~lambda333, 
                     mu1~0, q23~0, q32~0, q13~mu2, q12~mu3)

lnL3.0 <- make.classe(ttn3$tree, ttn3$states+1, 2)
lnL3.d <- make.classe(ttn3$tree, ttn3$states+1, 2, control=control.d)
lnL3.g <- make.classe(ttn3$tree, ttn3$states+1, 2, control=control.g)
lnL3.G <- make.classe(ttn3$tree, ttn3$states+1, 2, control=control.G)

lnL3.02 <- constrain(lnL3.0, lambda122~0, lambda211~0)
lnL3.d2 <- constrain(lnL3.d, lambda122~0, lambda211~0)
lnL3.g2 <- constrain(lnL3.g, lambda122~0, lambda211~0)
lnL3.G2 <- constrain(lnL3.G, lambda122~0, lambda211~0)

pars1.bisse <- c(1.4, 0.6, 0.2, 0.1, 0.5, 0.05)
pars1.classe <- c(rep(0, 6), pars1.bisse[-seq(2)])
names(pars1.classe) <- argnames(lnL1.0)
pars1.classe[c('lambda111', 'lambda222')] <- pars1.bisse[1:2]

pars2.geosse <- c(1.5, 0.5, 1.0, 0.7, 0.7, 1.4, 1.3)
names(pars2.geosse) <- argnames(lnL2.geosse)
pars2.classe <- 0 * starting.point.classe(ttn2$tree, 3)
pars2.classe['lambda222'] <- pars2.classe['lambda112'] <- pars2.geosse['sA']
pars2.classe['lambda333'] <- pars2.classe['lambda113'] <- pars2.geosse['sB']
pars2.classe['lambda123'] <-  pars2.geosse['sAB']
pars2.classe['mu2'] <- pars2.classe['q13'] <- pars2.geosse['xA']
pars2.classe['mu3'] <- pars2.classe['q12'] <- pars2.geosse['xB']
pars2.classe['q21'] <- pars2.geosse['dA']
pars2.classe['q31'] <- pars2.geosse['dB']
pars2.geosse2 <- pars2.geosse[c(3,1,2,4:7)]

pars3.classe <- c(5, 1, 0.1, 0.2, 2, 4, 3, 2.5, 2.1, 2.2)
names(pars3.classe) <- argnames(lnL3.0)
pars3.classe2 <- pars3.classe[-c(3,4)]

argvals <- list(condition.surv=TRUE)

f <- function(lik, pars, args)
  do.call(lik, c(list(pars), args))

ans <- -284.940610037178
expect_that(f(lnL1.b, pars1.bisse, argvals), equals(ans))
expect_that(f(lnL1.0, pars1.classe, argvals), equals(ans))
## expect_that(f(lnL1.d, pars1.classe, argvals), equals(ans))
## expect_that(f(lnL1.G, pars1.classe, argvals), equals(ans))
expect_that(f(lnL1.G, pars1.classe, argvals), equals(ans))
## 
expect_that(f(lnL1.02, pars1.bisse, argvals), equals(ans))
## expect_that(f(lnL1.d2, pars1.bisse, argvals), equals(ans))
## expect_that(f(lnL1.g2, pars1.bisse, argvals), equals(ans))
expect_that(f(lnL1.G2, pars1.bisse, argvals), equals(ans))

argvals <- list(condition.surv=FALSE)
ans <- -284.974405178804
expect_that(f(lnL1.b, pars1.bisse, argvals), equals(ans))
expect_that(f(lnL1.0, pars1.classe, argvals), equals(ans))
## expect_that(f(lnL1.d, pars1.classe, argvals), equals(ans))
## expect_that(f(lnL1.G, pars1.classe, argvals), equals(ans))
expect_that(f(lnL1.G, pars1.classe, argvals), equals(ans))
## 
expect_that(f(lnL1.02, pars1.bisse, argvals), equals(ans))
## expect_that(f(lnL1.d2, pars1.bisse, argvals), equals(ans))
## expect_that(f(lnL1.g2, pars1.bisse, argvals), equals(ans))
expect_that(f(lnL1.G2, pars1.bisse, argvals), equals(ans))

argvals <- list(condition.surv=TRUE, root=ROOT.GIVEN, root.p=c(0.6, 0.4))
ans <- -285.086707732767
expect_that(f(lnL1.b, pars1.bisse, argvals), equals(ans))
expect_that(f(lnL1.0, pars1.classe, argvals), equals(ans))
## expect_that(f(lnL1.d, pars1.classe, argvals), equals(ans))
## expect_that(f(lnL1.G, pars1.classe, argvals), equals(ans))
expect_that(f(lnL1.G, pars1.classe, argvals), equals(ans))
## 
expect_that(f(lnL1.02, pars1.bisse, argvals), equals(ans))
## expect_that(f(lnL1.d2, pars1.bisse, argvals), equals(ans))
## expect_that(f(lnL1.g2, pars1.bisse, argvals), equals(ans))
expect_that(f(lnL1.G2, pars1.bisse, argvals), equals(ans))

argvals <- list(condition.surv=TRUE, root=ROOT.EQUI)
ans <- -285.254147850837

expect_that(f(lnL1.b, pars1.bisse, argvals), equals(ans))
expect_that(f(lnL1.0, pars1.classe, argvals), equals(ans))
## expect_that(f(lnL1.d, pars1.classe, argvals), equals(ans))
## expect_that(f(lnL1.G, pars1.classe, argvals), equals(ans))
expect_that(f(lnL1.G, pars1.classe, argvals), equals(ans))
## 
expect_that(f(lnL1.02, pars1.bisse, argvals), equals(ans))
## expect_that(f(lnL1.d2, pars1.bisse, argvals), equals(ans))
## expect_that(f(lnL1.g2, pars1.bisse, argvals), equals(ans))
expect_that(f(lnL1.G2, pars1.bisse, argvals), equals(ans))

argvals <- list(condition.surv=TRUE)
ans <- -387.278039318393
expect_that(f(lnL2.geosse, pars2.geosse, argvals), equals(ans))
expect_that(f(lnL2.0, pars2.classe, argvals), equals(ans))
## expect_that(f(lnL2.d, pars2.classe, argvals), equals(ans))
## expect_that(f(lnL2.G, pars2.classe, argvals), equals(ans))
expect_that(f(lnL2.G, pars2.classe, argvals), equals(ans))
## 
expect_that(f(lnL2.02, pars2.geosse2, argvals), equals(ans))
## expect_that(f(lnL2.d2, pars2.geosse2, argvals), equals(ans))
## expect_that(f(lnL2.g2, pars2.geosse2, argvals), equals(ans))
expect_that(f(lnL2.G2, pars2.geosse2, argvals), equals(ans))

argvals <- list(condition.surv=FALSE, root=ROOT.GIVEN, root.p=c(0.5, 0.3, 0.2))
ans <- -386.963692579783
expect_that(f(lnL2.geosse, pars2.geosse, argvals), equals(ans))
expect_that(f(lnL2.0, pars2.classe, argvals), equals(ans))
## expect_that(f(lnL2.d, pars2.classe, argvals), equals(ans))
## expect_that(f(lnL2.G, pars2.classe, argvals), equals(ans))
expect_that(f(lnL2.G, pars2.classe, argvals), equals(ans))
## 
expect_that(f(lnL2.02, pars2.geosse2, argvals), equals(ans))
## expect_that(f(lnL2.d2, pars2.geosse2, argvals), equals(ans))
## expect_that(f(lnL2.g2, pars2.geosse2, argvals), equals(ans))
expect_that(f(lnL2.G2, pars2.geosse2, argvals), equals(ans))

argvals <- list(condition.surv=TRUE, root=ROOT.EQUI)
ans <- -387.179357960404
expect_that(f(lnL2.geosse, pars2.geosse, argvals), equals(ans))
expect_that(f(lnL2.0, pars2.classe, argvals), equals(ans))
## expect_that(f(lnL2.d, pars2.classe, argvals), equals(ans))
## expect_that(f(lnL2.G, pars2.classe, argvals), equals(ans))
expect_that(f(lnL2.G, pars2.classe, argvals), equals(ans))
## 
expect_that(f(lnL2.02, pars2.geosse2, argvals), equals(ans))
## expect_that(f(lnL2.d2, pars2.geosse2, argvals), equals(ans))
## expect_that(f(lnL2.g2, pars2.geosse2, argvals), equals(ans))
expect_that(f(lnL2.G2, pars2.geosse2, argvals), equals(ans))

argvals <- list(condition.surv=TRUE, root=ROOT.EQUI)
ans <- -387.179357960404
expect_that(f(lnL2.geosse, pars2.geosse, argvals), equals(ans))
expect_that(f(lnL2.0, pars2.classe, argvals), equals(ans))
## expect_that(f(lnL2.d, pars2.classe, argvals), equals(ans))
## expect_that(f(lnL2.G, pars2.classe, argvals), equals(ans))
expect_that(f(lnL2.G, pars2.classe, argvals), equals(ans))
## 
expect_that(f(lnL2.02, pars2.geosse2, argvals), equals(ans))
## expect_that(f(lnL2.d2, pars2.geosse2, argvals), equals(ans))
## expect_that(f(lnL2.g2, pars2.geosse2, argvals), equals(ans))
expect_that(f(lnL2.G2, pars2.geosse2, argvals), equals(ans))

ans <- 36.7578083566953
expect_that(lnL3.0(pars3.classe), equals(ans))
## expect_that(lnL3.d(pars3.classe), equals(ans))
## expect_that(lnL3.g(pars3.classe), equals(ans))
expect_that(lnL3.G(pars3.classe), equals(ans))

ans <- 35.5388000436388
expect_that(lnL3.02(pars3.classe2), equals(ans))
## expect_that(lnL3.d2(pars3.classe2), equals(ans))
## expect_that(lnL3.g2(pars3.classe2), equals(ans))
expect_that(lnL3.G2(pars3.classe2), equals(ans))

ans <- c(rep(3.05805572383494, 6),
         rep(4.58708358575242, 2),
         rep(4.58708358575242, 2) )
names(ans) <- argnames(lnL3.0)
expect_that(starting.point.classe(ttn3$tree, 2), equals(ans))

expect_that(diversitree:::stationary.freq.classe(pars1.classe, 2), 
            equals(diversitree:::stationary.freq.bisse(pars1.bisse)))

expect_that(diversitree:::stationary.freq.classe(pars2.classe, 3), 
            equals(diversitree:::stationary.freq.geosse(pars2.geosse)))

pars <- seq_len(27)  # for k = 3 states
names(pars) <- diversitree:::default.argnames.classe(3)
parlist <- diversitree:::inflate.pars.classe(pars, 3)
expect_that(parlist$lambda[3,1,2], equals(14))
expect_that(parlist$mu[2], equals(20))
expect_that(parlist$q[1,3], equals(23))
expect_that(diversitree:::flatten.pars.classe(parlist),
            is_identical_to(pars))

ans <- matrix(c(-78, 37, 43, 55, -87, 69, 81, 90, -99), nrow=3)
expect_that(diversitree:::projection.matrix.classe(pars, 3), equals(ans))

set.seed(1)
phy <- trees(pars3.classe, type="classe", n=2, max.t=1)
expect_that(lapply(phy, Ntip), equals(list(109, 8)))
lnL <- make.classe(phy[[1]], phy[[1]]$tip.state, 2)
expect_that(lnL(pars3.classe), equals(-7.29472293426202))

set.seed(3)
phy <- tree.classe(pars2.classe, max.t=3)
expect_that(as.numeric(table(phy$tip.state)), equals(c(19, 39, 17)))


