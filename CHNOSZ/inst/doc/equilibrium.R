## ----setup, include=FALSE, cache=FALSE---------------------------------------------
library(knitr)
## set global chunk options
opts_chunk$set(fig.align='center')
opts_chunk$set(tidy=TRUE)
## set code/output width
options(width=85)
## set number of digits
options(digits=3)

## ----AAsetup-----------------------------------------------------------------------
library(CHNOSZ)
data(thermo)
basis("CHNOS")
species(aminoacids(""))[1:5, ]

## ----AAfigures---------------------------------------------------------------------
res <- 200
aa <- aminoacids()

aaA <- function() {
  a <- affinity(O2=c(-90, -70, res), H2O=c(-20, 10, res))
  diagram(a, balance=1, names=aa)
}

aaB <- function() {
  a <- affinity(O2=c(-90, -70, res), H2O=c(-20, 10, res))
  e <- equilibrate(a, balance=1)
  diagram(e, names=aa)
}

aaC <- function() {
  a <- affinity(O2=c(-71, -66, res), H2O=c(-8, 4, res))
  diagram(a, balance="CO2", names=aa)
}

aaD <- function() {
  a <- affinity(O2=c(-71, -66), H2O=c(-8, 4))
  e <- equilibrate(a, balance="CO2")
  diagram(e, names=aa)
}

aaE <- function() {
  basis("O2", -66)
  a <- affinity(H2O=c(-8, 4))
  e <- equilibrate(a, balance="CO2")
  diagram(e, ylim=c(-5, -1), names=aa)
}

aaF <- function() {
  species(1:20, -4)
  a <- affinity(H2O=c(-8, 4))
  e <- equilibrate(a, balance="CO2")
  diagram(e, ylim=c(-5, -1), names=aa)
}

## ----AAnames-----------------------------------------------------------------------
AA <- aminoacids("")
names(AA) <- aa
AA

## ----showtime, include=FALSE-------------------------------------------------------
showtime <- function(st) {
  # plot time in lower-right of figure region
  f <- usrfig()
  par(xpd=TRUE)
  if(st[3] > 2) col <- "red" else col <- "black"
  text(f$x[2], f$y[1], paste(round(st[3], 1), "s\n"), adj=1, col=col)
  par(xpd=FALSE)
}

## ----AAplot, echo=FALSE, message=FALSE, results="hide", fig.width=13/2, fig.height=9/2, cache=TRUE----
layout(t(matrix(c(1:7, 11, 8:10, 12), nrow=4)), widths=c(1, 4, 4, 4), heights=c(1, 4, 4))

## row 0 (column titles)
opar <- par(mar=c(0, 0, 0, 0))
plot.new()
plot.new()
text(0.5, 0.5, "maximum affinity", cex=1.4)
plot.new()
text(0.5, 0.5, "equilibration", cex=1.4)
plot.new()
text(0.5, 0.5, "equilibration", cex=1.4)
par(opar)

## row 1 (balance = 1)
opar <- par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, "balance = 1", srt=90, cex=1.4)
par(opar)
# figure A
st <- system.time(dA <- aaA())
showtime(st)
title(main="loga(species) = -3", cex.main=1)
label.figure("A", col="blue", yfrac=0.9, xfrac=0.1)
# figure B
st <- system.time(dB <- aaB())
showtime(st)
title(main=paste("loga(total species) =", round(dB$loga.balance, 2)), cex.main=1)
label.figure("B", col="blue", yfrac=0.9, xfrac=0.1)

## row 2 (balance = nCO2)
opar <- par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, 'balance = "CO2"', srt=90, cex=1.4)
par(opar)
# figure C
st <- system.time(dC <- aaC())
showtime(st)
title(main="loga(species) = -3", cex.main=1)
label.figure("C", col="blue", yfrac=0.9, xfrac=0.1)
# figure D
st <- system.time(dD <- aaD())
showtime(st)
title(main=paste("loga(total CO2) =", round(dD$loga.balance, 2)), cex.main=1)
label.figure("D", col="blue", yfrac=0.9, xfrac=0.1)

## right (speciation at different total activity of CO2)
par(xpd=NA)
lines(c(-66, -64.5), c(4, 9), lty=2)
lines(c(-66, -64.5), c(-8, -8.5), lty=2)
par(xpd=FALSE)
# figure E
st <- system.time(dE <- aaE())
showtime(st)
title(main=paste("loga(total CO2) =", round(dE$loga.balance, 2)), cex.main=1)
label.figure("E", col="blue", yfrac=0.9, xfrac=0.1)
# figure F
st <- system.time(dF <- aaF())
showtime(st)
title(main=paste("loga(total CO2) =", round(dF$loga.balance, 2)), cex.main=1)
label.figure("F", col="blue", yfrac=0.9, xfrac=0.1)


## ----PRsetup, results="hide"-------------------------------------------------------
basis("CHNOS+")
organisms <- c("METJA", "HALJP", "METVO", "ACEKI", "GEOSE", "BACLI")
proteins <- c(rep("CSG", 3), rep("SLAP", 3))
species(proteins, organisms)

## ----PRfigures---------------------------------------------------------------------
prA <- function() {
  a <- affinity(O2=c(-90, -70, res), H2O=c(-20, 0, res))
  e <- equilibrate(a, balance="length", loga.balance=0)
  diagram(e, names=organisms)
}

prB <- function() {
  a <- affinity(O2=c(-90, -70))
  e <- equilibrate(a, balance="length", loga.balance=0)
  diagram(e, names=organisms, ylim=c(-5, -1), legend.x=NA)
}

prC <- function() {
  a <- affinity(O2=c(-90, -70, res), H2O=c(-20, 0, res))
  e <- equilibrate(a, normalize=TRUE, loga.balance=0)
  diagram(e, names=organisms)
}

prD <- function() {
  a <- affinity(O2=c(-90, -70))
  e <- equilibrate(a, normalize=TRUE, loga.balance=0)
  diagram(e, names=organisms, ylim=c(-5, -1), legend.x=NA)
}

prE <- function() {
  a <- affinity(O2=c(-90, -70, res), H2O=c(-20, 0, res))
  e <- equilibrate(a, as.residue=TRUE, loga.balance=0)
  diagram(e, names=organisms)
}

prF <- function() {
  a <- affinity(O2=c(-90, -70))
  e <- equilibrate(a, as.residue=TRUE, loga.balance=0)
  diagram(e, names=organisms, ylim=c(-3, 1), legend.x=NA)
}

## ----PRplot, message=FALSE, echo=FALSE, results="hide", fig.width=13/2, fig.height=9/2, cache=TRUE----
layout(t(matrix(1:12, nrow=4)), widths=c(1, 4, 4, 4), heights=c(1, 4, 4))

## row 0 (column titles)
opar <- par(mar=c(0, 0, 0, 0))
plot.new()
plot.new()
text(0.5, 0.5, 'balance = "length"', cex=1.4)
plot.new()
text(0.5, 0.5, "normalize = TRUE\n(balance = 1)", cex=1.4)
plot.new()
text(0.5, 0.5, "as.residue = TRUE\n(balance = 1)", cex=1.4)
par(opar)

## row 1 (equilibrate 2D)
opar <- par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, "equilibration", srt=90, cex=1.4)
par(opar)
# figure A (balance = "length")
st <- system.time(dA <- prA())
showtime(st)
label.figure("A", col="blue", yfrac=0.9, xfrac=0.1)
# figure C (normalize = TRUE)
st <- system.time(dC <- prC())
showtime(st)
label.figure("C", col="blue", yfrac=0.9, xfrac=0.1)
# figure E (as.residue = TRUE)
st <- system.time(dE <- prE())
showtime(st)
label.figure("E", col="blue", yfrac=0.9, xfrac=0.1)

## row 2 (equilibrate 1D)
opar <- par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, "equilibration", srt=90, cex=1.4)
par(opar)
# figure B (balance = "length")
st <- system.time(prB())
showtime(st)
label.figure("B", col="blue", yfrac=0.9, xfrac=0.1)
# figure D (normalize = TRUE)
st <- system.time(prD())
showtime(st)
label.figure("D", col="blue", yfrac=0.9, xfrac=0.1)
# figure F (as.residue = TRUE)
st <- system.time(prF())
showtime(st)
label.figure("F", col="blue", yfrac=0.9, xfrac=0.1)

## ----PRplot2, message=FALSE, results="hide", fig.width=6, fig.height=2, cache=TRUE----
layout(t(matrix(1:3)))
species(1:6, -111)
a <- affinity(O2=c(-90, -70, res), H2O=c(-20, 0, res))
d1 <- diagram(a, balance="length", names=organisms, main='balance = "length"')
d2 <- diagram(a, normalize=TRUE, names=organisms, main="normalize = TRUE")
d3 <- diagram(a, as.residue=TRUE, names=organisms, main="as.residue = TRUE")
stopifnot(!identical(d1$predominant, dA$predominant))
stopifnot(identical(d2$predominant, dC$predominant))
stopifnot(identical(d3$predominant, dE$predominant))

## ----Pourbaix, message=FALSE, results="hide", fig.width=4, fig.height=3, cache=TRUE----
basis(c("Cu+2", "H2O", "H+", "e-"))
species(c("Cu+2", "copper", "cuprite", "tenorite"))
for(loga in c(-1, 0, -2, -3)) {
  species("Cu+2", loga) 
  a <- affinity(pH=c(1.6, 7.6, 400), Eh=c(-0.2, 1, 400))
  if(loga==-1) d <- diagram(a)
  else d <- diagram(a, add=TRUE, names=NULL)
  iCu <- which(d$predominant == 1, arr.ind=TRUE)
  text(a$vals[[1]][max(iCu[, 1])] - 0.03, a$vals[[2]][min(iCu[, 2])] + 0.2, adj=1, loga)
}
water.lines()


## ----glutathione, results="hide", message=FALSE, fig.width=5, fig.height=3.5-------
basis(c("GSH", "NH3", "H2S", "H2O", "H+", "e-"))
basis("pH", 7)
species(c("GSH", "GSSG"))
a <- affinity(Eh=c(-0.3, -0.1))
# initial millimoles of GSH
mM <- c(10, 3, 1)
M <- mM * 1e-3
for(i in 1:3) {
  e <- equilibrate(a, loga.balance=log10(M[i]))
  diagram(e, alpha=TRUE, lty=c(0, i), add = i!=1, legend.x=NULL, ylim=c(0, 1), yline=1.6, lwd=2, ylab="aqueous species distribution")
  fGSH <- 1 - (10^e$loga.equil[[1]] / M[i])
  lines(e$vals[[1]], fGSH, col="blue", lty=i)
}
mtext(side=2, "fraction of GSH oxidized to GSSG", las=0, line=2.6, col="blue", cex=0.8)
mtext(side=2, "- - - - - - - - - - - - - - - - - - -", las=0, line=2.1, cex=0.8)
legend("topleft", lty=1:3, legend=paste(mM, "mM GSH"))

## ----Aksu_Doyle, message=FALSE, results="hide", fig.width=6, fig.height=4.8, cache=TRUE----
demo("copper", ask=FALSE)

## ----AddObigt----------------------------------------------------------------------
data(thermo)
add.obigt()

## ----ProteinFormation--------------------------------------------------------------
basis("CHNOS+")
species("CSG",c("METVO", "METJA"))

## ----ProteinInfo-------------------------------------------------------------------
protein.info(species()$name)

## ----ProteinAffinity---------------------------------------------------------------
a <- affinity()
a$values

## ----ProteinActivities-------------------------------------------------------------
e <- equilibrate(a)
e$loga.equil

## ----ProteinBasis------------------------------------------------------------------
protein.basis(species()$name, normalize=TRUE)

## ----ProteinEquil, results="hide"--------------------------------------------------
protein <- iprotein(c("CSG_METVO", "CSG_METJA"))
basis("CHNOS+")
swap.basis("O2", "H2")
protein.equil(protein, loga.protein=-3)

## ----ProteinSpeciation, results="hide", message=FALSE, fig.height=5.5, cache=TRUE----
organisms <- c("METSC", "METJA", "METFE", "HALJP",  "METVO", "METBU", "ACEKI", "GEOSE", "BACLI", "AERSA")
proteins <- c(rep("CSG", 6), rep("SLAP", 4))
basis("CHNOS+")
species(proteins, organisms)
a <- affinity(O2=c(-100, -65))
layout(matrix(1:2), heights=c(1, 2))
e <- equilibrate(a)
diagram(e, ylim=c(-2.8, -1.6), legend.x=NA, names=organisms)
water.lines(xaxis="O2")
title(main="Equilibrium activities of proteins, normalize = FALSE")
e <- equilibrate(a, normalize=TRUE)
diagram(e, ylim=c(-4, -2), legend.x=NA, names=organisms)
water.lines(xaxis="O2")
title(main="Equilibrium activities of proteins, normalize = TRUE")

## ----SulfurSpeciation, results="hide", message=FALSE, fig.width=6, fig.height=5.5----
basis("CHNOS+")
basis("pH", 5)
species(c("H2S", "S2-2", "S3-2", "S2O3-2", "S2O4-2",  "S3O6-2", "S5O6-2", "S2O6-2", "HSO3-", "SO2", "HSO4-"))
a <- affinity(O2=c(-50, -15), T=325, P=350)
layout(matrix(c(1, 2, 3, 3), nrow=2), widths=c(4, 1))
col <- rep(c("blue", "black", "red"), each=4)
lty <- 1:4
for(normalize in c(FALSE, TRUE)) {
  e <- equilibrate(a, loga.balance=-2, normalize=normalize)
  diagram(e, ylim=c(-30, 0), legend.x=NULL, col=col, lty=lty)
  water.lines(xaxis="O2", T=convert(325, "K"))
  title(main=paste("Aqueous sulfur speciation, normalize =", normalize))
}
par(mar=c(0, 0, 0, 0))
plot.new()
leg <- lapply(species()$name, expr.species)
legend("center", lty=lty, col=col, lwd=1.5, bty="n", legend=as.expression(leg), y.intersp=1.3)

