## ----setup, include=FALSE, cache=FALSE---------------------------------------------
library(knitr)
## set global chunk options
opts_chunk$set(fig.path='figure/hotspring-', cache.path='cache/hotspring-', fig.align='center', fig.show='hold', par=TRUE)
## set code/output width
options(width=85)
## set number of digits
options(digits=5)
## tune details of base graphics (http://yihui.name/knitr/hooks)
knit_hooks$set(par=function(before, options, envir){
if (before && options$fig.show!='none') par(mar=c(4,4,1,1),cex.lab=.95,cex.axis=.9,mgp=c(2,.7,0),tcl=-.3)
})

## ----libraryCHNOSZ-----------------------------------------------------------------
library(CHNOSZ)
data(thermo)
add.obigt()

## ----TpH---------------------------------------------------------------------------
bison.T <- c(93.3, 79.4, 67.5, 65.3, 57.1)
bison.pH <- c(7.350, 7.678, 7.933, 7.995, 8.257)

## ----TpHplot, fig.width=6, fig.height=3--------------------------------------------
distance <- c(0, 6, 11, 14, 22)
par(mfrow=c(1, 2), mar=c(4, 4, 3, 2))
xpoints <- seq(0, 22, length.out=128)
# T plot
plot(distance, bison.T, xlab="distance, m", ylab=axis.label("T"))
Tfun <- splinefun(distance, bison.T, method="mono")
lines(xpoints, Tfun(xpoints))

# pH plot
plot(distance, bison.pH, xlab="distance, m", ylab="pH")
pHfun <- splinefun(distance, bison.pH, method="mono")
lines(xpoints, pHfun(xpoints))

## ----proteins----------------------------------------------------------------------
# read the amino acid compositions
aa.annot <- read.aa(system.file("extdata/protein/DS11.csv", package="CHNOSZ"))
aa.phyla <- read.aa(system.file("extdata/protein/DS13.csv", package="CHNOSZ"))

## ----sitename----------------------------------------------------------------------
sites <- c("N", "S", "R", "Q", "P")
sitenames <- paste("bison", sites, sep="")

## ----classes-----------------------------------------------------------------------
classes <- unique(aa.annot$protein)
classes

## ----names-------------------------------------------------------------------------
# the names of the phyla in alphabetical order (except Deinococcus-Thermus at end)
phyla.abc <- sort(unique(aa.phyla$organism))[c(1:7,9:11,8)]
# an abbreviation for Dein.-Thermus
phyla.abbrv <- phyla.abc
phyla.abbrv[[11]] <- "Dein.-Thermus"
phyla.cols <- c("#f48ba5", "#f2692f", "#cfdd2a",
  "#962272", "#87c540", "#66c3a2", "#12a64a", "#f58656",
  "#ee3237", "#25b7d5", "#3953a4")
phyla.lty <- c(1:6, 1:5)
phyla.abbrv

## ----ZCplot, fig.width=5, fig.height=5, out.width='.49\\textwidth'-----------------
# 2011 plot
ylab <- expression(bar(italic(Z))[C])
plot(0, 0, xlim=c(-0.5, 5), ylim=c(-0.27, -0.11), xlab="location", xaxt="n", ylab=ylab)
axis(1, at=1:5)
col <- c("green", rep("black", 20))
lwd <- c(3, rep(1, 20))
clab <- c("hydrolase", "overall", "protease", "oxidoreductase", "transport", "membrane", "permease")
pf.annot <- protein.formula(aa.annot)
ZC.annot <- ZC(pf.annot)
for(i in 1:length(classes)) {
  lines(1:5, ZC.annot[(1:5)+5*(i-1)], col=col[i], lwd=lwd[i])
  if(classes[i] %in% clab) text(0.8, ZC.annot[1+5*(i-1)], classes[i], adj=1)
}
title(main="annotations")

# 2013 plot
pf.phyla <- protein.formula(aa.phyla)
ZC.phyla <- ZC(pf.phyla)
# set up plot
plot(0, 0, xlim=c(1, 5), ylim=c(-0.27, -0.11), xlab="location", ylab=ylab)
for(i in 1:length(phyla.abc)) {
  # which of the model proteins correspond to this phylum
  iphy <- which(aa.phyla$organism==phyla.abc[i])
  # the locations (of 1, 2, 3, 4, 5) where this phylum is found
  ilocs <- match(aa.phyla$protein[iphy], sitenames)
  # the plotting symbol: determined by alphabetical position of the phylum
  points(ilocs, ZC.phyla[iphy], pch=i-1, cex=1.2)
  # a line to connect same phyla occurring at adjacent sites
  inlocs <- rep(NA, 5)
  inlocs[ilocs] <- ilocs
  lines(inlocs, ZC.phyla[iphy][match(1:5, ilocs)])
}
legend("bottomright", pch=0:10, legend=phyla.abbrv, bg="white", cex=0.9)
title(main="major phyla")

## ----setup.basis-------------------------------------------------------------------
setup.basis <- function() {
  basis(c("HCO3-", "H2O", "NH3", "HS-", "H2", "H+"))
  basis(c("HCO3-", "NH3", "HS-", "H+"), c(-3, -4, -7, -7.933))
}

## ----overall-----------------------------------------------------------------------
setup.basis()
ip.annot <- add.protein(aa.annot)
species("overall", sitenames)

## ----residue-----------------------------------------------------------------------
pl <- protein.length(ip.annot[1:5])
mysp <- species()
mysp[, 1:6]/pl

## ----logaH2------------------------------------------------------------------------
get.logaH2 <- function(T) -11 + T * 3/40

## ----affinity----------------------------------------------------------------------
species(1:5, 0)
a <- affinity(T=bison.T, pH=bison.pH, H2=get.logaH2(bison.T))
a.res <- t(as.data.frame(a$values))/pl
a.res

## ----affinitymax-------------------------------------------------------------------
apply(a.res, 2, which.max)

## ----affinityadj-------------------------------------------------------------------
a.res <- a.res - log10(pl)
apply(a.res, 2, which.max)

## ----Tlim--------------------------------------------------------------------------
Tlim <- c(50, 100)

## ----TlogaH2plot, fig.width=6, fig.height=3, message=FALSE, results='hide'---------
par(mfrow=c(1, 2))
# first plot
a <- affinity(T=Tlim, H2=c(-7, -4))
diagram(a, fill=NULL, names=1:5, normalize=TRUE)
lines(Tlim, get.logaH2(Tlim), lty=3)
# second plot
species(1:5, -3)
xT <- Tfun(xpoints)
xpH <- pHfun(xpoints)
xH2 <- get.logaH2(xT)
a <- affinity(T=xT, pH=xpH, H2=xH2)
a$vars[1] <- "distance, m"
a$vals[[1]] <- xpoints
e <- equilibrate(a, normalize=TRUE)
diagram(e, legend.x=NULL)
legend("bottom", lty=1:5, legend=1:5, bty="n", cex=0.6)

## ----stripplot, fig.width=6, fig.height=3.5, message=FALSE, results='hide'---------
loadclass <- function(class) {
  species(delete=TRUE)
  species(rep(class, each=5), rep(sitenames, length(class)))
}
xclasses <- c("overall", "transferase", "transport", "synthetase", "membrane", "permease")
loadclass(xclasses)
a <- affinity(T=xT, pH=xpH, H2=xH2)
a$vars[1] <- "distance, m"
a$vals[[1]] <- xpoints
col <- c("red", "orange", "yellow", "green", "blue")
par(mfrow=c(1, 2), mar=c(4, 4, 1, 1))
for(i in 1:2) {
  ispecies <- lapply((1:3)+(i-1)*3, function(x) {1:5+(x-1)*5} )
  names(ispecies) <- xclasses[(1:3)+(i-1)*3]
  strip(a = a, ispecies = ispecies, col = col, xticks = distance, cex.names = 1)
}

## ----methionine, fig.width=6, fig.height=4, message=FALSE, results='hide', cache=TRUE----
par(mfrow=c(2, 3))
for(j in 1:2) {
  # use old [Met] for first row and new [Met] for second row
  if(j==2) {
    data(thermo)
    ip.annot <- add.protein(aa.annot)
  }
  # setup basis species and proteins
  setup.basis()
  # make the plots
  for(annot in c("overall", "transferase", "synthase")) {
    ip <- ip.annot[aa.annot$protein==annot]
    a <- affinity(T=c(50, 100), H2=c(-7, -4), iprotein=ip)
    diagram(a, fill=NULL, names=1:5, normalize=TRUE)
    # add logaH2-T line
    lines(par("usr")[1:2], get.logaH2(par("usr")[1:2]), lty=3)
    # add a title
    title(main=annot)
  }
}

## ----alpha.blast-------------------------------------------------------------------
alpha.blast <- function() {
  out <- xtabs(ref ~ protein + organism, aa.phyla)
  # put it in correct order, then turn counts into fractions
  out <- out[c(1,5:2), c(1:7,9:11,8)]
  out <- out/rowSums(out)
  return(out)
}

## ----alpha.equil-------------------------------------------------------------------
alpha.equil <- function(i=1) {
  # order the names and counts to go with the alphabetical phylum list
  iloc <- which(aa.phyla$protein==sitenames[i])
  iloc <- iloc[order(match(aa.phyla$organism[iloc], phyla.abc))]
  # set up basis species, with pH specific for this location
  setup.basis()
  basis("pH", bison.pH[i])
  # calculate metastable equilibrium activities of the residues
  a <- affinity(H2=c(-11, -1, 101), T=bison.T[i], iprotein=ip.phyla[iloc])
  e <- equilibrate(a, loga.balance=0, as.residue=TRUE)
  # remove the logarithms to get relative abundances
  a.residue <- 10^sapply(e$loga.equil, c)
  colnames(a.residue) <- aa.phyla$organism[iloc]
  # the BLAST profile
  a.blast <- alpha.blast()
  # calculate Gibbs energy of transformation (DGtr) and find optimal logaH2
  iblast <- match(colnames(a.residue), colnames(a.blast))
  r <- revisit(e, "DGtr", log10(a.blast[i, iblast]), plot.it=FALSE)
  # return the calculated activities, logaH2 range, DGtr values, and optimal logaH2
  return(list(alpha=a.residue, H2vals=a$vals[[1]], DGtr=r$H, logaH2.opt=r$xopt))
}

## ----alphaplot, fig.width=8, fig.height=6, tidy=FALSE, message=FALSE, results='hide', cache=TRUE----
ip.phyla <- add.protein(aa.phyla)
layout(matrix(1:6, ncol=3), heights=c(2, 1))
equil.results <- list()
for(i in 1:5) {
  # get the equilibrium degrees of formation and the optimal logaH2
  ae <- alpha.equil(i)
  equil.results[[i]] <- ae
  if(i %in% c(1, 3, 5)) {
    iphy <- match(colnames(ae$alpha), phyla.abc)
    # top row: equilibrium degrees of formation
    thermo.plot.new(xlim=range(ae$H2vals), ylim=c(0, 0.5), xlab=axis.label("H2"),
      ylab=expression(alpha[equil]), yline=2, cex.axis=1, mgp=c(1.8, 0.3, 0))
    for(j in 1:ncol(ae$alpha)) {
      lines(ae$H2vals, ae$alpha[, j], lty=phyla.lty[iphy[j]])
      ix <- seq(1, length(ae$H2vals), length.out=11)
      ix <- head(tail(ix, -1), -1)
      points(ae$H2vals[ix], ae$alpha[, j][ix], pch=iphy[j]-1)
    } 
    title(main=paste("site", i))
    legend("topleft", pch=iphy-1, lty=phyla.lty[iphy], legend=phyla.abbrv[iphy], bg="white")
    # bottom row: Gibbs energy of transformation and position of minimum
    thermo.plot.new(xlim=range(ae$H2vals), ylim=c(0, 1/log(10)), xlab=axis.label("H2"),
      ylab=expr.property("DGtr/2.303RT"), yline=2, cex.axis=1, mgp=c(1.8, 0.3, 0))
    lines(ae$H2vals, ae$DGtr)
    abline(v=ae$logaH2.opt, lty=2)
    abline(v=get.logaH2(bison.T[i]), lty=3, lwd=1.5)
    if(i==1) legend("bottomleft", lty=c(3, 2), lwd=c(1.5, 1), bg="white",
      legend=c("Equation 2", "optimal"))
  }
}

## ----E.AgAgCl----------------------------------------------------------------------
E.AgAgCl <- function(T) {
  0.23737 - 5.3783e-4 * T - 2.3728e-6 * T^2 - 2.2671e-9 * (T+273)
}

## ----meters------------------------------------------------------------------------
T.ORP <- c(93.9, 87.7, 75.7, 70.1, 66.4, 66.2)
pH.ORP <- c(8.28, 8.31, 7.82, 7.96, 8.76, 8.06)
ORP <- c(-258, -227, -55, -58, -98, -41)

## ----ORP2Eh, message=FALSE, results='hide'-----------------------------------------
Eh <- ORP/1000 + E.AgAgCl(T.ORP)
pe <- convert(Eh, "pe", T=convert(T.ORP, "K"))
logK.ORP <- subcrt(c("e-", "H+", "H2"), c(-2, -2, 1), T=T.ORP)$out$logK
logaH2.ORP <- logK.ORP - 2*pe - 2*pH.ORP

## ----sulfur, message=FALSE, results='hide'-----------------------------------------
loga.HS <- log10(c(4.77e-6, 2.03e-6, 3.12e-7, 4.68e-7, 2.18e-7))
loga.SO4 <- log10(c(2.10e-4, 2.03e-4, 1.98e-4, 2.01e-4, 1.89e-4))
logK.S <- subcrt(c("HS-", "H2O", "SO4-2", "H+", "H2"), c(-1, -4, 1, 1, 4), T=bison.T)$out$logK
logaH2.S <- (logK.S + bison.pH - loga.SO4 + loga.HS) / 4

## ----oxygen, message=FALSE, results='hide'-----------------------------------------
DO <- c(0.173, 0.776, 0.9, 1.6, 2.8)
logaO2 <- log10(DO/1000/32)
logK <- subcrt(c("O2", "H2", "H2O"), c(-0.5, -1, 1), T=bison.T)$out$logK
logaH2.O <- 0 - 0.5*logaO2 - logK

## ----logaH2plot, fig.width=4, fig.height=4, out.width='.45\\textwidth', cache=TRUE----
# 2011 plot
xlab <- axis.label("T")
ylab <- axis.label("H2")
plot(Tlim, get.logaH2(Tlim), xlim=Tlim, ylim=c(-45,0),
  xlab=xlab, ylab=ylab, type="l", lty=3)
points(T.ORP, logaH2.ORP, pch=15)
lines(T.ORP, logaH2.ORP, lty=2)
points(bison.T, logaH2.O, pch=16)
lines(bison.T, logaH2.O, lty=2)
points(bison.T, logaH2.S, pch=17)
lines(bison.T, logaH2.S, lty=2)
llab <- c("Equation 2", "ORP", "dissolved oxygen", "sulfate/sulfide")
text(c(65, 80, 80, 74), c(-4, -25, -40, -11), llab)

# 2013 plot
plot(Tlim, get.logaH2(Tlim), xlim=Tlim, ylim=c(-11,-2),
  xlab=xlab, ylab=ylab, type="l", lty=3)
lines(bison.T, sapply(equil.results, "[", "logaH2.opt"), lty=2)
points(bison.T, sapply(equil.results, "[", "logaH2.opt"), pch=21, bg="white")
text(90, -5.3, "Equation 2")
text(66, -9, "optimal parameterization\nfor observed\nphylum abundances", adj=0)

## ----alphaplot5, fig.width=9, fig.height=6, out.width='.88\\textwidth', cache=TRUE----
layout(matrix(c(1, 2, 3, 4, 5, 6), nrow=2, byrow=TRUE), widths=c(2, 2, 2))
par(mar=c(2.5, 0, 2.5, 0))
plot.new()
legend("topright", pch=0:11, legend=phyla.abbrv, bty="n", cex=1.5)
lim <- c(-6, -0.5)
equil.opt <- a.blast <- alpha.blast()
for(iloc in 1:5) {
  a.equil <- equil.results[[iloc]]
  iopt <- match(a.equil$logaH2.opt, a.equil$H2vals)
  ae.opt <- a.equil$alpha[iopt, ]
  # which are these phyla in the alphabetical list of phyla
  iphy <- match(names(ae.opt), phyla.abc)
  equil.opt[iloc, iphy] <- ae.opt
  mar <- c(2.5, 4.0, 2.5, 1)
  thermo.plot.new(xlab=expression(log[2]*alpha[obs]), ylab=expression(log[2]*alpha[equil]),
    xlim=lim, ylim=lim, mar=mar, cex=1, yline=1.5)
  # add points and 1:1 line
  points(log2(a.blast[iloc, iphy]), log2(ae.opt), pch=iphy-1)
  lines(lim, lim, lty=2)
  title(main=paste("site", iloc))
  # within-plot legend: DGtr
  DGexpr <- as.expression(quote(Delta*italic(G[tr])/italic(RT) == phantom()))
  DGval <- format(round(2.303*a.equil$DGtr[iopt], 3), nsmall=3)
  legend("bottomright", bty="n", legend=c(DGexpr, DGval))
}

## ----barchart, fig.width=6, fig.height=3-------------------------------------------
par(mar=c(4, 4, 3, 0), mgp=c(1.8, 0.7, 0))
par(mfrow=c(1, 3), cex=1)
# make the blast plot
ab <- alpha.blast()
rownames(ab) <- 1:5
barplot(t(ab), col=phyla.cols, ylab=NULL, xlab="site", axes=TRUE, cex.axis=0.8, cex.names=0.8, las=1)
mtext(expression(alpha[obs]), 2, 2, cex=1.1*par("cex"))
title(main="BLAST profile", cex.main=0.8)
# make the equilibrium plot
rownames(equil.opt) <- 1:5
barplot(t(equil.opt), col=phyla.cols, ylab=NULL, xlab="site", axes=TRUE, cex.axis=0.8,
  cex.names=0.8, las=1)
mtext(expression(alpha[equil]), 2, 2, cex=1.1*par("cex"))
title(main="metastable\nequilibrium", cex.main=0.8)
# add legend
par(mar=c(4, 1, 3, 0))
plot.new()
legend("bottomleft", legend=rev(phyla.abbrv), fill=rev(phyla.cols), bty="n", cex=0.7)

