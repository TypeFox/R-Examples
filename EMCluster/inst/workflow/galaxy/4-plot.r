rm(list = ls(all.names = TRUE))

library(EMCluster)
load("./data/galaxy.em.rda")

x.hist <- hist(x, nclass = 100, plot = FALSE)

ret.plot <- NULL
for(K in 1:6){
  if(ret[[K]]$em$llhdval > ret[[K]]$Rnd$llhdval){
    ret.plot[[K]] <- ret[[K]]$em
  } else{
    ret.plot[[K]] <- ret[[K]]$Rnd
  }  
}

xlim <- c(8000, 36000)
xx <- seq(xlim[1], xlim[2], length = 200)
kmax <- 6

yy <- list()
for(K in 1:kmax){
  yy[[K]] <- rep(0, length(xx))
  for(k in 1:K){
    yy[[K]] <- yy[[K]] + ret.plot[[K]]$pi[k] *
               dnorm(xx, ret.plot[[K]]$Mu[k, 1],
                         sqrt(ret.plot[[K]]$LTSigma[k, 1]))
  }
}
ylim <- range(do.call("c", lapply(yy, range)))

postscript("plot/galaxy.ps", width = 6, height = 6, horizontal = FALSE)
par(mar = c(5, 4, 4, 4), mgp = c(2, 1, 0))
hist(x, nclass = 50, xlim = xlim, main = "Corona Borealis",
     xlab = "Velocities of 82 Galaxies", ylab = "Frequency")
box()
for(K in kmax:1){
  lines(x = xx, y = yy[[K]] * 11 / (ylim[2] - ylim[1]),
        lty = K, col = K, lwd = 2)
}
axis(4, at = seq(0, 11, length = 5),
     labels = formatC(seq(0, ylim[2], length = 5), format = "e", digits = 1))
mtext("Density", side = 4, line = 2)
legend(28000, 7, paste("K=", 1:kmax, sep = ""),
       lty = 1:kmax, col = 1:kmax, lwd = 2)
dev.off()
