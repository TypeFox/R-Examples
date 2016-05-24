# Fe-minerals and aqueous species in Fe-S-O-H-C system
# after Garrels and Christ, 1965 Figure 7.21
# to reproduce their diagram as closely as posssible, use their thermodynamic data (from Appendix 2)
mod.obigt(c("Fe+2", "Fe+3"), G=c(-20300, -2520))
mod.obigt(c("hematite", "magnetite", "pyrrhotite", "pyrite", "siderite"), G=c(-177100, -242400, -23320, -36000, -161060))
mod.obigt(c("SO4-2", "HS-", "H2S", "HSO4-"), G=c(-177340, 3010, -6540, -179940))
mod.obigt(c("CO2", "HCO3-", "CO3-2"), G=c(-92310, -140310, -126220))
# conditions and system definition
pH <- c(0, 14, 400)
Eh <- c(-1, 1, 400)
T <- 25
basis(c("FeO", "SO4-2", "H2O", "H+", "e-", "CO3-2"))
basis("SO4-2", -6)
basis("CO3-2", 0)
species(c("Fe+2", "Fe+3"), -6)
species(c("pyrrhotite", "pyrite", "hematite", "magnetite", "siderite"))
# two sets of changing basis species:
# speciate SO4-2, HSO4-, HS-, H2S as a function of Eh and pH
# speciate CO3-2, HCO3-, CO2 as a function of pH
bases <- c("SO4-2", "HSO4-", "HS-", "H2S")
bases2 <- c("CO3-2", "HCO3-", "CO2")
# calculate affinities using the predominant basis species
# using blend=TRUE we get curvy lines, particularly at the boundaries with siderite
# compare with the plot in Garrels and Christ, 1965
m1 <- mosaic(bases, bases2, TRUE, pH=pH, Eh=Eh, T=T)
# make a diagram and add water stability lines
diagram(m1$A.species)
water.lines("pH", "Eh", T=convert(T, "K"), col="seagreen", lwd=1.5)
# show lines for Fe(aq) = 10^-4 M
species(c("Fe+2", "Fe+3"), -4)
m2 <- mosaic(bases, bases2, TRUE, pH=pH, Eh=Eh, T=T)
diagram(m2$A.species, add=TRUE, names=NULL, dotted=3)
title(main=paste("Iron oxides, sulfides and carbonate in water, log(total S) = -6,",
  "log(total C)=0, after Garrels and Christ, 1965", sep="\n"))
# overlay the carbonate basis species predominance fields
diagram(m1$A.bases2, add=TRUE, col="blue", col.names="blue", dotted=3)
# reset the database, as it was changed in this example
data(thermo)
