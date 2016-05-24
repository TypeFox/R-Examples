### R code from vignette source 'RDieHarder.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library(RDieHarder)
options(SweaveHooks=list(twofig=function() {par(mfrow=c(1,2))},
                         twofig2=function() {par(mfrow=c(2,1))},
                         onefig=function() {par(mfrow=c(1,1))}))
prettyVersion <- packageDescription("RDieHarder")$Version
prettyDate <- format(Sys.Date(), "%B %e, %Y")


###################################################
### code chunk number 2: loaddata
###################################################
if (file.exists("RDieHarder.Rdata")) load("RDieHarder.Rdata")


###################################################
### code chunk number 3: rd-example
###################################################
  if (!exists("dh")) dh <- dieharder("ran0","diehard_2dsphere",seed=2)
  #dh
  plot(dh)


###################################################
### code chunk number 4: rd-example1
###################################################
  if (!exists("dh1")) dh1 <- dieharder("mt19937","diehard_2dsphere",seed=2)
  #dh1
  plot(dh1)


###################################################
### code chunk number 5: r-rngs
###################################################
rngs <- c("R_wichmann_hill", "R_marsaglia_multic",
          "R_super_duper", "R_mersenne_twister",
          "R_knuth_taocp", "R_knuth_taocp2")

if (!exists("rl")) rl <- lapply(rngs, function(rng) dieharder(rng, "diehard_runs", seed=12345))

oldpar <- par(mfrow=c(2,3), mar=c(2,3,3,1))
invisible(lapply(rl, function(res) {
  qqplot(res$data, seq(0, 1, length.out=length(res$data)),
         main=paste(res$generator, ":", round(res$p.value, digits=3)),
         ylab="", type="S")
  abline(0, 1, col='gray')
}))
par(oldpar) # reset graph defaults



###################################################
### code chunk number 6: <savedata
###################################################
save(dh, dh1, rl, file="RDieHarder.Rdata")


