#!/usr/bin/env r

suppressMessages(library(RDieHarder))

rngs <- c("R_wichmann_hill", "R_marsaglia_multic", 
          "R_super_duper", "R_mersenne_twister",
          "R_knuth_taocp", "R_knuth_taocp2")

rl <- lapply(rngs, function(rng) dieharder(rng, "runs", seed=12345))

oldpar <- par(mfrow=c(2,3), mar=c(2,3,3,1))
invisible( lapply(rl, function(res) {
  qqplot(res$data, seq(0, 1, length.out=length(res$data)),
         main=paste(res$generator, ":", round(res$p.value, digits=4)),
         ylab="", type="S")
  abline(0, 1, col='gray')
}))
par(oldpar)
