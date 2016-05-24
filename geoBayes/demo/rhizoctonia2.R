######################################################################
## 
### Commentary: EB analysis of the rhizoctonia data. Uses the 
### functions mcsglmm, bf1skel, bf2new, bf2optim.
## 
######################################################################

library(geoBayes)


### Load data
data(rhizoctonia)


### Create prediction grid
predgrid <- mkpredgrid2d(rhizoctonia[c("Xcoord", "Ycoord")],
                         par.x = 100, chull = TRUE, exf = 1.2)
rhizdata <- stackdata(rhizoctonia, predgrid$grid)


### Define the model
corrf <- "spherical"
kappa <- 0
ssqdf <- 1
ssqsc <- 1
betm0 <- 0
betQ0 <- .01
linkp <- "probit"


### Skeleton points
philist <- c(100, 140, 180)
omglist <- c(0, .5, 1, 1.5)
parlist <- expand.grid(phi=philist, linkp=linkp, omg=omglist, kappa = kappa)
estimate <- list(linkp = linkp, phi = c(100, 200), omg = c(0, 2),
                 kappa = kappa)


### MCMC sizes
Nout <- 1000
Nthin <- 10
Nbi <- 300


### Take MCMC samples
runs <- list()
for (i in 1:NROW(parlist)) {
  runs[[i]] <- mcsglmm(Infected ~ 1, 'binomial', rhizdata, weights = Total,
                       atsample = ~ Xcoord + Ycoord,
                       Nout = Nout, Nthin = Nthin, Nbi = Nbi,
                       betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
                       phistart = parlist$phi[i], omgstart = parlist$omg[i],
                       linkp = parlist$linkp[i], kappa = parlist$kappa[i], 
                       corrfcn = corrf, phisc = 0, omgsc = 0)
}

bf <- bf1skel(runs)

bfall <- bf2new(bf, phi = seq(100, 200, 10), omg = seq(0, 2, .2))


plotbf2(bfall, c("phi", "omg"))

plotbf2(bfall, c("phi", "omg"), profile = TRUE, type = "b", ylab="log(BF)")

bf2optim(bf, estimate)

