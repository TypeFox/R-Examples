## Buffer + ionization: Metastablilities of
## thiol peroxidases from model bactera
## (ECOLI, BACSU mesophile; AQUAE thermophile,
## THIDA acidophile, BACHD alkaliphile)
basis("CHNOS+")
organisms <- c("ECOLI", "AQUAE", "BACSU", "BACHD", "THIDA")
species("TPX", organisms)
# create a buffer with our proteins in it
mod.buffer("TPX", paste("TPX", organisms, sep="_"))
# set up the buffered activities
basis(c("CO2", "H2O", "NH3", "O2"), "TPX")
a <- affinity(return.buffer=TRUE, T=50)
basis(c("CO2", "H2O", "NH3", "O2"), as.numeric(a[1:4]))
a <- affinity(pH=c(4, 10, 300), T=c(40, 60, 300))
e <- equilibrate(a, normalize=TRUE)
diagram(e, fill=NULL)
title(main="Thiol peroxidases from bacteria")
legend("topleft", describe.basis(thermo$basis[-6,]))

## Buffer + ionization: relative stabilities
## of E. coli sigma factors on a T-pH diagram
# (sigma factors 24, 32, 38, 54, 70, i.e.
# RpoE, RpoH, RpoS, RpoN, RpoD)
proteins <- c("RPOE", "RP32", "RPOS", "RP54", "RPOD")
basis("CHNOS+")
basis("pH", 7.4)
# define and set the buffer
mod.buffer("sigma", paste(proteins, "ECOLI", sep="_"))
basis(c("CO2", "NH3", "H2S", "O2"), "sigma")
logact <- affinity(return.buffer=TRUE, T=25)
# Set the activities of the basis species to constants 
# corresponding to the buffer, and diagram the relative
# stabilities as a function of T and pH
basis(c("CO2", "NH3", "H2S", "O2"), as.numeric(logact))
species(paste(proteins, "ECOLI", sep="_"))
a <- affinity(pH=c(5, 10), T=c(10, 40))
diagram(a, normalize=FALSE, fill="heat")
title(main="Relative stabilities of sigma factors in E. coli")
ptext <- c(describe.property("T", 25), 
  describe.basis(ibasis=c(2, 6), oneline=TRUE))
btext <- describe.basis(ibasis=c(1, 3, 4, 5), oneline=TRUE)
legend("bottomleft", legend=c("preset (input values):",
  ptext, "buffered (results):", btext), bty="n")
