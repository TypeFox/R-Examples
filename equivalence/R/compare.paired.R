"compare.paired" <-
function(x, y=NULL, techniques=
c("tost", "ptte", "twos", "sign", "srank", "mwhit", "boot", "manly"),
Epsilon=1, reps=50)
{
### VALIDATE INPUT
   cat("Techniques:\n")
   print(techniques)
   if(is.character(y)) {
      cat("CAUTION. In place of y, techniques taken as:", deparse(y), "\n")
      techniques <- y
      y <- NULL
   }

   if(is.null(y)) {
      d <- x
      ntechniques <- length(techniques)
      techniques <- techniques[techniques != "twos" 
                             & techniques != "mwhit"]
      if(ntechniques != length(techniques)) {
         cat("Single sample received.  Techniques:\n"); print(techniques)
      }
   } else d <- y-x

   x <- x[!is.na(x)]
   y <- y[!is.na(y)]
   d <- d[!is.na(d)]
   n <- length(d)

   if("sign" %in% techniques) {
      p1 <- p2 <- 0.5
      repeat{
         p1 <- p1 - 0.01
         p2 <- 1 - p1
         mid <- as.integer(n/2)
         if(pbinom(mid, n, p2) < 0.0001 | p2 > 0.97) break
      }
   }


  cat("Estimating powers by bootstrap ...\n")
  par(las=1)
  if(n < 50) Epsilon <- 1 else Epsilon <- 0.5
  Epsilon.tost <- Epsilon * sd(d)
  npoints <- 15
  mu <- -2.0*Epsilon.tost + 4*Epsilon.tost * (0:npoints)/npoints
  power.twos <- power.ptte <- power.tost <- rep(NA, npoints+1)
  power.sign <- power.srank <- power.mwhit <- power.tost
  power.boot <- power.manly <- power.tost

  for(i in 1:(npoints+1)) {

   if("tost" %in% techniques) {
    power.tost[i] <- mean(boot(d-mean(d)+mu[i], tost.boot, reps,
                            Epsilon=Epsilon.tost)$t)
   }

   if("ptte" %in% techniques) {
    power.ptte[i] <- mean(boot(d-mean(d, na.rm=TRUE)+mu[i], ptte.boot, reps,
                            Epsilon=Epsilon)$t)
   }

   if("sign" %in% techniques) {
    power.sign[i] <- mean(boot(d-mean(d, na.rm=TRUE)+mu[i], sign.boot, reps,
                            p1=p1, p2=p2)$t)
   }

   if("srank" %in% techniques) {
    power.srank[i] <- mean(boot(d-mean(d, na.rm=TRUE)+mu[i], srank.boot, reps)$t)
   }

   if("twos" %in% techniques) {
    xy <- matrix(c(x-mean(x), y-mean(y)+mu[i]), ncol=2)
    power.twos[i] <- mean(boot(xy, twos.boot, reps,
                            Epsilon=Epsilon*sqrt(2))$t)
   }

   if("mwhit" %in% techniques) { "x & y INTERCHANGED: d = y-x"
    yx <- matrix(c(y-mean(y)+mu[i], x-mean(x)), ncol=2)
    power.mwhit[i] <- mean(boot(yx, mwhit.boot, reps)$t)
   }

   if("boot" %in% techniques) {
    power.boot[i] <- mean(boot(d-mean(d)+mu[i], boot.boot, reps,
                            Epsilon=Epsilon.tost)$t)
    }

   if("manly" %in% techniques) {
    power.manly[i] <- mean(boot(d-mean(d)+mu[i], manly.boot, reps,
                            Epsilon=Epsilon.tost)$t)
    }
  }

##### PLOT POWERS ########################################################
  maxim <- 1
  plot( mu, power.tost, axes=FALSE, type="n",
        xlab=expression(paste(mu[D])),
        ylab="Power", 
        xlim=2*c(-Epsilon.tost, Epsilon.tost), ylim=c(0, 1.1) )
  box()
  ntypes <- length(techniques)
  lines(mu, power.tost,   col="blue", lty = 1, lwd=3)
  lines(mu, power.ptte,    col="red", lty = 2, lwd=3)
  lines(mu, power.twos, col="purple", lty = 3, lwd=3)
  lines(mu, power.sign,  col="brown", lty = 4, lwd=3)
  lines(mu, power.srank,col="orange", lty = 5, lwd=3)
  lines(mu, power.mwhit, col="green", lty = 6, lwd=3)
  lines(mu, power.boot,  col="black", lty = 3, lwd=3)
  lines(mu, power.manly,  col="pink", lty = 2, lwd=3)
  abline(h=0.05, lty=2); abline(h=maxim, lty=2)
  axis(2, c(0.05, maxim))
  axis(1, c(-Epsilon.tost, Epsilon.tost), c("-Epsilon", paste("Epsilon =",
  as.numeric(formatC(Epsilon.tost, format = "f", digits = 2)) )))
  abline(v=-Epsilon.tost, lty=2); abline(v=Epsilon.tost, lty=2)

### Order: "tost", "ptte", "twos", "sign", "srank", "mwhit", "boot", "manly"
  colours <- NULL
  ltys <- NULL
  if("tost" %in% techniques) { colours <- c(colours, "blue")
                               ltys <- c(ltys, 1) }
  if("ptte" %in% techniques) { colours <- c(colours, "red")
                               ltys <- c(ltys, 2) }
  if("twos" %in% techniques) { colours <- c(colours, "purple")
                               ltys <- c(ltys, 3) }
  if("sign" %in% techniques) { colours <- c(colours, "brown")
                               ltys <- c(ltys, 4) }
  if("srank" %in% techniques) { colours <- c(colours, "orange")
                               ltys <- c(ltys, 5) }
  if("mwhit" %in% techniques) { colours <- c(colours, "green")
                               ltys <- c(ltys, 6) }
  if("boot" %in% techniques) { colours <- c(colours, "black")
                               ltys <- c(ltys, 3) }
  if("manly" %in% techniques) { colours <- c(colours, "pink")
                               ltys <- c(ltys, 2) }
  legend(Epsilon.tost, 1.0, techniques, col=colours, lwd=rep(2, ntypes),
                                lty=ltys)
#############################################################################

  cat("n =", n, "Sign test equivalence region:\n")
  print(c( as.numeric(formatC(p1, format = "f", digits = 2)),
                      as.numeric(formatC(p2, format = "f", digits = 2)) ))

  cat("Signed rank equivalence region: pnorm( +/- 0.5*sqrt(2))\n")
  print(c( as.numeric(formatC(pnorm(-0.5*sqrt(2)), format = "f", digits = 2)), 
                            as.numeric(formatC(pnorm(+0.5*sqrt(2)),
                               format = "f", digits = 2)) ))

  cat("Mann-Whitney equivalence region: pnorm(-.5/sqrt(2))~pnorm(1/sqrt(2))\n")
  print(c( as.numeric(formatC(pnorm(-0.5/sqrt(2)), format = "f", digits = 2)), 
                            as.numeric(formatC(pnorm(+1.0/sqrt(2)),
                               format = "f", digits = 2)) ))

  data.frame(MU=mu, TOST=power.tost, PTTE=power.ptte, TWOS=power.twos, 
             SIGN=power.sign, SIGNRANK=power.srank, MWHIT=power.mwhit,
             MANLY=power.manly)
}

