######################################################################
## 
### Commentary: EB analysis of the rhizoctonia data using the GEV
### link. 
## 
######################################################################

library(geoBayes)

data(rhizoctonia)

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
linkp <- seq(-2, 2, .5)

### Skeleton points
philist <- c(100, 150, 200)
omglist <- 2
parlist <- expand.grid(phi=philist, linkp=linkp, omg=omglist, kappa = kappa)
estimate <- list(linkp = c(-2, 2), phi = philist, omg = omglist,
                 kappa = kappa)

### MCMC sizes
Nout <- Npro <- 100
Nthin <- Nprt <- 10
Nbi <- Nprb <- 30


### Take MCMC samples
runs <- list()
for (i in 1:NROW(parlist)) {
  runs[[i]] <- mcsglmm(Infected ~ 1, 'GEV.binomial', rhizdata, weights = Total,
                       atsample = ~ Xcoord + Ycoord,
                       Nout = Nout, Nthin = Nthin, Nbi = Nbi,
                       betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
                       phistart = parlist$phi[i], omgstart = parlist$omg[i],
                       linkp = parlist$linkp[i], kappa = parlist$kappa[i], 
                       corrfcn = corrf, phisc = 0, omgsc = 0)
}
bf <- bf1skel(runs, transf = TRUE)
bfall <- bf2new(bf, phi = seq(50, 250, 10))
plotbf2(bfall, c("phi", "linkp"), 1)

est <- ebsglmm(Infected ~ 1, 'GEV.binomial', rhizdata, weights = Total,
               atsample = ~ Xcoord + Ycoord, parskel = parlist,
               paroptim = estimate, corrfcn = corrf, 
               Nout = Nout, Nthin = Nthin, Nbi = Nbi,
               Npro = Npro, Nprt = Nprt, Nprb = Nprb, 
               betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
               dispersion = 1, useCV=TRUE, transf = TRUE)

est$parest
