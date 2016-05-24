### R code from vignette source 'anintro.Rnw'

###################################################
### code chunk number 1: anintro.Rnw:25-26
###################################################
  if(exists(".orig.enc")) options(encoding = .orig.enc)


###################################################
### code chunk number 2: anintro.Rnw:55-56
###################################################
options(width=90)


###################################################
### code chunk number 3: install_CHNOSZ (eval = FALSE)
###################################################
## install.packages("CHNOSZ")


###################################################
### code chunk number 4: library_CHNOSZ
###################################################
library(CHNOSZ)


###################################################
### code chunk number 5: data_thermo
###################################################
data(thermo)


###################################################
### code chunk number 6: info_ethylene
###################################################
info("ethylene")


###################################################
### code chunk number 7: info_88
###################################################
info(88)


###################################################
### code chunk number 8: info_ethylene_gas
###################################################
info("ethylene","gas")


###################################################
### code chunk number 9: info_acetic
###################################################
aadata <- info(info("acetic acid"))
print(aadata)


###################################################
### code chunk number 10: summary_thermo
###################################################
summary(thermo)


###################################################
### code chunk number 11: browse.refs (eval = FALSE)
###################################################
## browse.refs(88)


###################################################
### code chunk number 12: info_acid
###################################################
info("acid")


###################################################
### code chunk number 13: info_spaceacid
###################################################
info(" acid")


###################################################
### code chunk number 14: info_OH
###################################################
info("(OH)")


###################################################
### code chunk number 15: protein_LYSC
###################################################
ip <- iprotein("LYSC_CHICK")
aa <- ip2aa(ip)
aa2eos(aa)


###################################################
### code chunk number 16: formula_LYSC
###################################################
pf <- protein.formula(aa)
as.chemical.formula(pf)


###################################################
### code chunk number 17: info_LYSC
###################################################
si <- info("LYSC_CHICK")
info(si)


###################################################
### code chunk number 18: subcrt_water
###################################################
subcrt("water")


###################################################
### code chunk number 19: subcrt_C2H5OH
###################################################
subcrt(c("C2H5OH","O2","CO2","H2O"),c(-1,-3,2,3),T=37)


###################################################
### code chunk number 20: subcrt_C2H5OH_O2aq
###################################################
subcrt(c("C2H5OH","O2","CO2","H2O"),c(-1,-3,2,3),c("aq","aq","aq","liq"),T=37)


###################################################
### code chunk number 21: subcrt_C2H5OH_unbal
###################################################
subcrt(c("C2H5OH","CO2","H2O"),c(-1,2,3),T=37)


###################################################
### code chunk number 22: basis_not (eval = FALSE)
###################################################
## basis(c("CO2","H2O","NH3","H2S","H+"))


###################################################
### code chunk number 23: basis
###################################################
basis(c("CO2","H2O","NH3","O2","H2S","H+"))


###################################################
### code chunk number 24: subcrt_C2H5OH_auto
###################################################
subcrt(c("C2H5OH","CO2","H2O"),c(-1,2,3),T=37)


###################################################
### code chunk number 25: subcrt_C2H5OH
###################################################
subcrt(c("C2H5OH"),c(-1),T=37)


###################################################
### code chunk number 26: subcrt_ethanol
###################################################
subcrt(c("ethanol","acetaldehyde"),c(-1,1),T=37)


###################################################
### code chunk number 27: subcrt_acetaldehyde
###################################################
data(thermo)
basis(c("glutamic acid","methionine","isoleucine","lysine","tyrosine","H+"))
subcrt(c("ethanol","acetaldehyde"),c(-1,1),T=37)


###################################################
### code chunk number 28: subcrt_LYSC
###################################################
data(thermo)
basis("CHNOS+")
subcrt("LYSC_CHICK",1,T=25)


###################################################
### code chunk number 29: subcrt_ALAT1 (eval = FALSE)
###################################################
## aa <- uniprot.aa("ALAT1_HUMAN")
## add.protein(aa)
## subcrt("ALAT1_HUMAN",1,T=25)


###################################################
### code chunk number 30: Bjerrum_diagram
###################################################
basis("CHNOS+")
species(c("CO2", "HCO3-", "CO3-2"))
a <- affinity(pH=c(4, 12))
e <- equilibrate(a)
diagram(e, ylim=c(-6, 0))
a <- affinity(pH=c(4, 12), T=150)
e <- equilibrate(a)
diagram(e, add=TRUE, col="red")


###################################################
### code chunk number 31: demo_solubility
###################################################
demo("solubility", ask=FALSE)


###################################################
### code chunk number 32: CSG_species
###################################################
basis("CHNOS")
species(c("SLAP_ACEKI", "CSG_METJA", "CSG_METVO", "CSG_HALJP"))


###################################################
### code chunk number 33: CSG_affinity
###################################################
a <- affinity(O2=c(-90,-70))


###################################################
### code chunk number 34: CSG_diagram
###################################################
e <- equilibrate(a, normalize=TRUE)
diagram(e, legend.x="bottomleft", ylim=c(-6, -2))


###################################################
### code chunk number 35: PredominanceDiagram
###################################################
species(c("SLAP_ACEKI", "SLAP_GEOSE", "SLAP_BACLI", "SLAP_AERSA"))
basis(c("NH3", "H2S"), c(-1, -10))
a <- affinity(O2=c(-85, -70), H2O=c(-5, 0))
diagram(a, normalize=TRUE)


###################################################
### code chunk number 36: residue.info
###################################################
protein.basis(species()$name, normalize=TRUE)


###################################################
### code chunk number 37: Bowers
###################################################
basis(c("HCl","H2O","Ca+2","CO2","Mg+2","SiO2","O2","H+"),
  c(999,0,999,999,999,999,999,-7))
species(c("quartz","talc","forsterite","tremolite","diopside",
  "wollastonite","monticellite","merwinite"))
a <- affinity("Mg+2"=c(-12,-4),"Ca+2"=c(-8,0),T=300,P=1000)
diagram(a)


###################################################
### code chunk number 38: example_diagram (eval = FALSE)
###################################################
## example(diagram)


###################################################
### code chunk number 39: examples (eval = FALSE)
###################################################
## examples()


###################################################
### code chunk number 40: demo (eval = FALSE)
###################################################
## demo("findit")


###################################################
### code chunk number 41: helpthermo (eval = FALSE)
###################################################
## help(thermo)


###################################################
### code chunk number 42: amylaseplot
###################################################
basis("CHNOSe")
basis(c("NH3", "H2S"), c(-6, -3))
species(c("AMY_BACSU", "O08452_PYRFU"))
a <- affinity(pH=c(2, 10), Eh=c(-1, 1))
diagram(a, balance="CO2", lwd=2, fill=NULL, cex.names=1.5)
a <- affinity(pH=c(2, 10), Eh=c(-1, 1), T=100)
diagram(a, balance="CO2", lwd=2, fill=NULL, col="red", names=NULL, add=TRUE)
water.lines()
water.lines(T=398.15, col="red")
BKMdat <- read.csv(system.file("extdata/cpetc/BKM60_Fig7.csv", package="CHNOSZ"))
points(BKMdat$pH, BKMdat$Eh, pch=20)
points(c(8.5, 6.2, 4.2, 6.7), c(0.018, 0.223, 0.022, 0.067), pch=3, col="red")
points(c(8.24, 7.17, 8.11), c(-0.44, -0.41, -0.49), pch=17, col="red")
ltext <- c(describe.property("T", 25), describe.property("T", 100))
legend("topright", legend=ltext, lty=1, col=c("black", "red"))
ltext <- c("soils [BKM60]", "yellowstone [SWMP05]", "iceland [SA02]")
legend("bottomleft", legend=ltext, pch=c(20, 3, 17))


###################################################
### code chunk number 43: yeastplot
###################################################
locations <- yeastgfp()
gfp <- yeastgfp(locations)
aa <- more.aa(gfp$protein, "Sce")
for(i in 1:length(locations)) {
  avgaa <- aasum(aa[[i]], gfp$abundance[[i]], average=TRUE, protein=locations[i])
  add.protein(avgaa)
}
basis("CHNOS+")
species(locations, "Sce")
a <- affinity(O2=c(-82, -65))
e <- equilibrate(a, loga.balance=0, normalize=TRUE)
mycolor <- topo.colors(length(locations))
diagram(e, names=locations, ylim=c(-5, -3), legend.x=NA,
  col=mycolor, lwd=2)
dp <- describe.property(c("T", "P"), c(25, 1))
db <- describe.basis(ibasis=(1:6)[-5])
legend("topright", legend=c(dp, db), bty="n")


###################################################
### code chunk number 44: bufferplot
###################################################
layout(matrix(1:2, nrow=1), widths=c(2, 1))
b.species <- c("Fe", "CO2", "H2O", "N2", "H2", "H2S", "SiO2")
b.state <- c("cr1", "gas", "liq", "gas", "gas", "aq", "aq")
b.logact <- c(0, 1, 0, 0, 0, 0, 0)
basis(b.species, b.state, b.logact)
xlim <- c(0, 350)
thermo.plot.new(xlim=xlim, ylim=c(-4, 4), xlab=axis.label("T"), ylab=axis.label("H2"))
bufferline <- function(buffer, ixlab) {
  basis("H2", buffer)
  a <- affinity(T=xlim, P=300, return.buffer=TRUE, exceed.Ttr=TRUE)
  lines(a$vals[[1]], a$H2)
  text(a$vals[[1]][ixlab], a$H2[ixlab], buffer)
}
bufferline("FeFeO", 20)
bufferline("QFM", 38)
bufferline("PPM", 102)
bufferline("HM", 51)
basis("H2", 0)
for(logact in c(-6, -10, -15)) {
  species(c("formaldehyde", "HCN"), logact)
  a <- affinity(T=xlim, P=300)
  d <- diagram(a, what="H2", lty=c(2, 3), add=TRUE)
  text(a$vals[[1]][13], mean(sapply(d$plotvals, c)[13, ]), logact)
}
plot.new()
legend("topleft", legend = c(describe.property("P", 300), describe.basis(ibasis=c(2,4)),
  "minerals", "HCN", "formaldehyde"), lty=c(NA, NA, NA, 1, 2, 3), bg="white")


###################################################
### code chunk number 45: revisit
###################################################
basis("CHNOS")
species(c("isoleucine", "tyrosine", "glutamic acid", "methionine", "aspartic acid"))
a <- affinity(CO2=c(2, 5), O2=c(-72, -70))
e <- equilibrate(a, balance=1)
r <- revisit(e)
title(main=paste("CV minimum =", round(r$optimum, 2)))


###################################################
### code chunk number 46: revisit_alpha
###################################################
basis(c("CO2", "O2"), c(r$x, r$y))
a <- affinity()
par(mfrow=c(1, 2))
e <- equilibrate(a, balance=1)
d <- diagram(e, alpha=TRUE, names=aminoacids(1, species()$name))
plot.new()
legend("topleft", describe.basis(basis()), bg="white")


###################################################
### code chunk number 47: findit
###################################################
basis("CHNOS")
species(c("isoleucine", "tyrosine", "glutamic acid", "methionine", "aspartic acid"))
f <- findit(list(CO2=c(-5, 5), O2=c(-85, -65), H2S=c(-10, 5), H2O=c(-10, 0)), 
  niter=5, res=10, balance=1)


###################################################
### code chunk number 48: findit_alpha
###################################################
a <- affinity()
par(mfrow=c(1, 2))
e <- equilibrate(a, balance=1)
d <- diagram(e, alpha=TRUE, names=aminoacids(1, species()$name))
plot.new()
legend("topleft", describe.basis(basis()), bg="white")


###################################################
### code chunk number 49: session_info
###################################################
sessionInfo()


