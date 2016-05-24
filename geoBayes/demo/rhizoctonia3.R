######################################################################
## 
### Commentary: EB analysis of the rhizoctonia data.
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
linkp <- "probit"

### Skeleton points
philist <- c(100, 140, 180)
omglist <- c(0, .5, 1, 1.5)
parlist <- expand.grid(phi=philist, linkp=linkp, omg=omglist, kappa = kappa)
estimate <- list(linkp = linkp, phi = c(100, 200), omg = c(0, 2),
                 kappa = kappa)

### MCMC sizes
Nout <- Npro <- 1000
Nthin <- Nprt <- 10
Nbi <- Nprb <- 300

est <- ebsglmm(Infected ~ 1, 'binomial', rhizdata, weights = Total,
               atsample = ~ Xcoord + Ycoord, parskel = parlist,
               paroptim = estimate, corrfcn = corrf, 
               Nout = Nout, Nthin = Nthin, Nbi = Nbi,
               Npro = Npro, Nprt = Nprt, Nprb = Nprb, 
               betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
               dispersion = 1, useCV=TRUE)

library(geoR)
z0pred <- rowMeans(est$mcmcsample$z[!est$mcmcsample$whichobs, ])
geoR:::image.kriging(locations = expand.grid(predgrid$xycoord),
                     borders = predgrid$borders, values = z0pred,
                     x.leg = c(3150, 3450), y.leg = c(1120, 1200),
                     col = gray((64:32)/64))
