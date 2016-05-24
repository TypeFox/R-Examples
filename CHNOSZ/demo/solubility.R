# demo showing how to calculate CO2(gas) or calcite solubility and aqueous carbonate speciation
# the affinity() ... equilibrate() sequence in CHNOSZ gives *metastable equilibrium activities*
#   (activities for a total activity of the balanced component given in loga.balance) ...
# here we are interested in finding the value of loga.balance itself
#   (the total activity of [CO2, HCO3-, CO3-2] species in the aqueous phase).
# this total activity is the solubility of CO2(gas) or calcite if the affinities of the
#   aqueous species as formed from CO2(gas) or calcite are all equal to zero.
# note that the affinities for species in metastable equilibrium are all equal.
#   Afun() calculates the metastable equilibrium affinities for a given loga.balance
#   and uniroot() finds the loga.balance where they are zero
# additionally, if we are reacting calcite, the activity of Ca+2 should be set equal to loga.balance

# for comparison with published calcite solubility plot, see Fig. 4A in
# Manning et al., 2013, Reviews in Mineralogy & Geochemistry, v. 75, pp. 109-148
# (doi: 10.2138/rmg.2013.75.5)

# for comparison with published CO2 solubility plot, see Fig. 4.5 in
# Stumm and Morgan, 1996, Aquatic Chemistry: Chemical Equilibria and Rates in Natural Waters
# (New York: John Wiley & Sons), 3rd edition

# set this to CO2 or calcite
#what <- "CO2"
what <- "calcite"

# function to return the affinity of the metastable equilibrium species
Afun <- function(loga.balance=-3, T=25) {
  if(what=="calcite") basis("Ca+2", loga.balance)
  a <- affinity(T=T)
  e <- equilibrate(a, loga.balance=loga.balance)
  # set metastable activities and re-calculate the affinity
  species(1:3, unlist(e$loga.equil))
  a <- affinity(T=T)
  # check they're actually equal
  stopifnot(all(abs(unlist(a$values) - a$values[[1]]) < 1e-10))
  return(a$values[[1]])
}

# set up system
if(what=="CO2") {
  basis("CHNOS+")
  basis("CO2", "gas")
  # ca. atmospheric PCO2
  basis("CO2", -3.5)
} else if(what=="calcite") {
  basis(c("calcite", "Ca+2", "H2O", "O2", "H+"))
}
species(c("CO2", "HCO3-", "CO3-2"))
T <- 25
# decrease this for higher resolution
pHstep <- 1

# where we'll store the results
loga.tot <- numeric()
loga.CO2 <- loga.HCO3 <- loga.CO3 <- numeric()

# loop over pH range
pHs <- seq(0, 14, pHstep)
for(pH in pHs) {
  print(paste("pH =", pH))
  basis("pH", pH)
  # this is for the solubility
  loga.balance <- suppressMessages(uniroot(Afun, c(-10, 10), T=T)$root)
  loga.tot <- c(loga.tot, loga.balance)
  # this is for the speciation
  if(what=="calcite") basis("Ca+2", loga.balance)
  a <- affinity(T=T)
  e <- equilibrate(a, loga.balance=loga.balance)
  loga <- unlist(e$loga.equil)
  loga.CO2 <- c(loga.CO2, loga[1])
  loga.HCO3 <- c(loga.HCO3, loga[2])
  loga.CO3 <- c(loga.CO3, loga[3])
}

# make plot
ylim <- c(-10, 4)
thermo.plot.new(xlim=range(pHs), ylim=ylim, xlab="pH", ylab="log a")
lines(pHs, loga.tot, lwd=4, col="green2")
lines(pHs, loga.CO2, lwd=2)
lines(pHs, loga.HCO3, lty=2, lwd=2)
lines(pHs, loga.CO3, lty=3, lwd=2)
legend(ifelse(what=="calcite", "topright", "topleft"), lty=c(1, 1:3), lwd=c(4, 2, 2, 2), col=c("green2", rep("black", 3)),
       legend=as.expression(c("total", expr.species("CO2", state="aq"), expr.species("HCO3-"), expr.species("CO3-2"))))
title(main=substitute(what~"solubility at"~T~degree*"C", list(what=what, T=T)))
