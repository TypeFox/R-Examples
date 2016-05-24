## NaCl dissocation logK f(T,P)
## after Shock et al., 1992, Fig. 1
## (Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A. and Helgeson, H. C. (1992) 
##  Calculation of the thermodynamic properties of aqueous species at high pressures and temperatures: 
##  Effective electrostatic radii, dissociation constants and standard partial molal properties to 1000 degrees C and 5 kbar. 
##  J. Chem. Soc. Faraday Trans. 88, 803-826. http://dx.doi.org/10.1039/FT9928800803 )
species <- c("NaCl", "Na+", "Cl-")
coeffs <- c(-1, 1, 1)
# start a new plot and show the experimental logK
thermo.plot.new(xlim=c(0, 1000), ylim=c(-5.5, 1),
  xlab=axis.label("T"), ylab=axis.label("logK"))
expt <- read.csv(system.file("extdata/cpetc/SOJSH.csv", 
  package="CHNOSZ"), as.is=TRUE)
points(expt$T,expt$logK, pch=expt$pch)
# we'll be at 9 distinct pressure conditions, including Psat
P <- c(list("Psat"), as.list(seq(500, 4000, by=500)))
# for each of those what's the range of temperature
T <- list()
# T > 350 degC at Psat is possibly inappropriate; see "Warning" of subcrt.Rd 
T[[1]] <- seq(0, 370, 5)
T[[2]] <- seq(265, 465, 5)
T[[3]] <- seq(285, 760, 5)
T[[4]] <- seq(395, 920, 5)
T[[5]] <- T[[6]] <- T[[7]] <- T[[8]] <- T[[9]] <- seq(400, 1000, 5)
# calculate and plot the logK
logK <- numeric()
for(i in 1:length(T)) {
  s <- subcrt(species, coeffs, T=T[[i]], P=P[[i]])
  lines(s$out$T, s$out$logK)
  # keep the calculated values for each experimental condition
  iexpt <- which(P[[i]]==expt$P)
  Texpt <- expt$T[iexpt]
  logK <- c(logK, splinefun(s$out$T, s$out$logK)(Texpt))
}
legend("bottomleft",pch=unique(expt$pch),
  legend=c(unique(expt$source),tail(expt$source,1)))
title(main=paste("NaCl(aq) = Na+ + Cl-\n",
  "Psat and 500-4000 bar, after Shock et al., 1992"))
# where do we diverge most from experiment?
imaxdiff <- which.max(abs(logK - expt$logK))
stopifnot(all.equal(c("Psat", 347.7),
  as.character(expt[imaxdiff,1:2])))
# what's our average divergence?
stopifnot(mean(abs(logK - expt$logK)) < 0.09)
