######################################################################
## 
### Commentary: Estimate the effect of the rhizoctonia disease on
### yield using EB methods.
## 
######################################################################

library(geoBayes)

### Load the data
data(rhizoctonia)
rhiz <- na.omit(rhizoctonia)
rhiz$IR <- rhiz$Infected/rhiz$Total # Incidence rate of the
                              # rhizoctonia disease

### Define the model
corrf <- "spherical"
ssqdf <- 1
ssqsc <- 1
tsqdf <- 1
tsqsc <- 1
betm0 <- 0
betQ0 <- diag(.01, 2, 2)

### Skeleton points
philist <- seq(120, 280, 40)
linkp <- 1
omglist <- seq(0, 2.5, .5)
parlist <- expand.grid(phi = philist, linkp = linkp, omg = omglist)
estimate <- list(linkp = linkp, phi = c(100, 300), omg = c(0, 2))


### MCMC sizes
Nout <- Npro <- 1000
Nthin <- Nprt <- 10
Nbi <- Nprb <- 3000


est <- ebstrga(Yield ~ IR, rhiz, 
               atsample = ~ Xcoord + Ycoord, parskel = parlist,
               paroptim = estimate, corrfcn = corrf, 
               Nout = Nout, Nthin = Nthin, Nbi = Nbi,
               Npro = Npro, Nprt = Nprt, Nprb = Nprb, 
               betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
               tsqdf = tsqdf, tsqsc = tsqsc,
               useCV = TRUE)

est$parest

mceb <- mcmcmake(est$mcmcsample)

library(mcmcplots)
ipar <- grep(paste(c("phi", "tsq", "ssq", "beta(_[0-9]+)?", "omg"),
                   collapse = "|"), dimnames(mceb)[[2]])

traplot(mceb[, ipar])
denplot(mceb[, ipar])
summary(mceb[, ipar])


