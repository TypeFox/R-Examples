opar <- par(ask = interactive() &&
            (.Device %in% c("X11", "GTK", "gnome", "windows","quartz")))

data(phyto)

# smooth the observed data
bobs  <- dpill(obs$t, obs$y)
n     <- tail(obs$t, n = 1) - obs$t[1] + 1
obss  <- ksmooth(obs$t, obs$y, kernel = "normal", bandwidth = bobs, 
                 n.points = n)
obss  <- as.data.frame(obss)
names(obss) <- c("t", "y")
obss  <- na.omit(obss[match(sim$t, obss$t),])


# plot data
basedate <- as.Date("1960/1/1")
plot(basedate + sim$t, sim$y, type="l", xaxt="n",
  xlim=c(basedate + min(obs$t, sim$t), basedate + max(obs$t, sim$t)),
  ylim=c(min(obs$y, sim$y), max(obs$y, sim$y)),
  xlab="1994", ylab="Phytoplankton (mg/L)",
  col="grey60", font=2, lwd=2, cex.lab=1.4, las=1)
lines(basedate + obss$t, obss$y, lwd=2)
points(basedate + obs$t, obs$y, lwd=2)

# positions and names to draw x-axes
lab.mon  <- as.Date(c(paste("1994/", c(1,3,5,7,9,11), "/01", sep=""), "1995/1/1"))
lab.lab  <- c("Jan","Mar","May","Jul","Sep","Nov","Jan")
axis(1, at=lab.mon, labels=lab.lab, font=2, cex.lab=1.4)

legend(as.Date("1994/09/01"), 50, c("measurement", "smoothed measurement", "simulation"), 
lty=c(0,1,1), pch=c(1,NA,NA), lwd=2, col=c(1,1,"grey60"))


# semiquantitative comparison

sqc <- compareME(obs$y, sim$y, obs$t, sim$t, type = c("normalized"))
round(sqc$normalized, digits=3)
compareME(type="name")


# time transformation

tt <- timeTransME(obss$y, sim$y, obss$t, sim$t, ME = SMSE, time = "transform", 
    type = "normalized", trials = 5, timeME = MAE, timeMEtype = "dissimilarity")
cat("Error:", round(tt$totalME, digits = 3),"\n")
cat("Time stress:", round(tt$timeME, digits = 3),"\n")

# plot time transformation

basedate <- as.Date("1960/1/1")
plot(basedate + sim$t, sim$y, type="l", xaxt="n",
  ylim=c(min(obs$y, sim$y), max(obs$y, sim$y)), 
  xlab="time", ylab="Phytoplankton (mg/L)",
  col="grey60", font=2, lwd=2, cex.lab=1.2, las=1)
lines(basedate + obss$t, obss$y, lwd=2)
points(basedate + obs$t, obs$y, lwd=2)
lines(basedate + tt$x, tt$yp, lwd=2, col="grey60", lty=2)

# positions and names to draw x-axes
lab.mon  <- as.Date(c(paste("1994/", c(1,3,5,7,9,11), "/01", sep=""), "1995/1/1"))
lab.lab  <- c("Jan","Mar","May","Jul","Sep","Nov","Jan")
axis(1, at=lab.mon, labels=lab.lab, font=2, cex.lab=1.4)

legend(basedate + 12600, 50, c("measurement", "smoothed measurement", "simulation",
"transformed simulation"), lty=c(0,1,1,2), pch=c(1,NA,NA,NA), lwd=2, col=c(1,1,"grey60","grey60"))


# interval sequences

obs_f1 <- f.slope(obss$t, obss$y)
sim_f1 <- f.slope(sim$t,  sim$y)
lcs    <- LCS(obs_f1, sim_f1)
cat("QSI:",lcs$QSI,"\n")

qvlcs <- qvalLCS(obs$y, sim$y, obs$t, sim$t, smooth = "obs")
plot.qvalLCS(qvlcs)

par(opar)
