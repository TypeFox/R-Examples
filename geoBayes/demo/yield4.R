######################################################################
## 
### Commentary: Estimate the effect of the rhizoctonia disease on
### yield using EB methods. Also estimates the link function.
## 
######################################################################

library(geoBayes)

### Load the data
data(rhizoctonia)
rhiz <- na.omit(rhizoctonia)
rhiz$IR <- rhiz$Infected/rhiz$Total # Incidence rate of the
                              # rhizoctonia disease
scale <- exp(mean(log(rhiz$Yield)))
rhiz$YieldScaled <- rhiz$Yield/scale

### Define the model
corrf <- "spherical"
ssqdf <- 1
ssqsc <- 1/20
tsqdf <- 1
tsqsc <- 1
betm0 <- 0
betQ0 <- diag(.01, 2, 2)

### Skeleton points
philist <- seq(120, 280, 40)
linkp <- c(0, .3, .5, .7, 1)
omglist <- seq(0, 1, .25)
parlist <- expand.grid(phi = philist, linkp = linkp, omg = omglist)
estimate <- list(linkp = c(0, 1), phi = c(100, 300), omg = c(0, 1))


### MCMC sizes
Nout <- Npro <- 1000
Nthin <- Nprt <- 10
Nbi <- Nprb <- 300


### Estimation
est <- ebstrga(YieldScaled ~ IR, rhiz, 
               atsample = ~ Xcoord + Ycoord, parskel = parlist,
               paroptim = estimate, corrfcn = corrf, 
               Nout = Nout, Nthin = Nthin, Nbi = Nbi,
               Npro = Npro, Nprt = Nprt, Nprb = Nprb, 
               betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
               tsqdf = tsqdf, tsqsc = tsqsc,
               useCV = TRUE)

### Parameter estimates
est$parest

mceb <- mcmcmake(est$mcmcsample)

ipar <- grep(paste(c("phi", "tsq", "ssq", "beta(_[0-9]+)?", "omg"),
                   collapse = "|"), dimnames(mceb)[[2]])

library(mcmcplots)
traplot(mceb[, ipar])
denplot(mceb[, ipar])
summary(mceb[, ipar])

