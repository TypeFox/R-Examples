######################################################################
## 
### Commentary: MCMC analysis of the rhizoctonia data.
## 
######################################################################

library(geoBayes)

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
phiprior <- c(100, 1, 1000, 100) # U(100, 200)
phisc <- 4
omgprior <- c(2, 1, 1, 0)        # Exp(mean = 2)
omgsc <- .3
linkp <- "probit"

### MCMC sizes
Nout <- 1000
Nthin <- 10
Nbi <- 300

emt <- mcsglmm(Infected ~ 1, 'binomial', rhizdata, weights = Total,
               atsample = ~ Xcoord + Ycoord,
               Nout = Nout, Nthin = Nthin, Nbi = Nbi,
               betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
               phipars = phiprior, omgpars = omgprior, linkp = linkp, 
               corrfcn = corrf, kappa = kappa, phisc = phisc, omgsc = omgsc, 
               dispersion = 1, test = TRUE)

emc <- mcsglmm(Infected ~ 1, 'binomial', rhizdata, weights = Total,
               atsample = ~ Xcoord + Ycoord,
               Nout = Nout, Nthin = Nthin, Nbi = Nbi,
               betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
               phipars = phiprior, omgpars = omgprior, linkp = linkp, 
               corrfcn = corrf, kappa = kappa, phisc = phisc, omgsc = omgsc, 
               dispersion = 1, test = FALSE)

emcmc <- mcmcmake(emc)

plot.ts(cbind(phi = emc$phi, omg = emc$omg, beta = c(emc$beta), ssq = emc$ssq),
        nc = 2)

summary(emcmc[, c("phi", "omg", "beta", "ssq")])

plot(emcmc[, c("phi", "omg", "beta", "ssq")])

library(geoR)
z0pred <- rowMeans(emc$z[!emc$whichobs, ])
geoR:::image.kriging(locations = predgrid$xygrid,
                     borders = predgrid$borders, values = z0pred,
                     x.leg = c(3150, 3450), y.leg = c(1120, 1200),
                     col = gray((64:32)/64))

