## Nucleobase - Amino Acid Interaction Eh-H2O
# for this example we try a custom basis set
basis(c("CO2", "H2O", "glutamine", "e-", "H+"), c(-3, 0, -3, 0, -7))
species(c("uracil", "cytosine", "adenine", "guanine",
  "phenylalanine", "proline", "lysine", "glycine"), "aq")
# this loaded four nucleobases and four related amino acids
# (coded for by the homocodon triplets)
# check out the predominance diagrams
a.1 <- affinity(H2O=c(-5, 0, 300), Eh=c(-0.5, 0, 300))
diagram(a.1, fill=NULL)
# overlay a different temperature
a.2 <- affinity(H2O=c(-5, 0, 300), Eh=c(-0.5, 0, 300), T=100)
diagram(a.2, col="red", add=TRUE, names=NULL)
# add title and legend
title(main="Nucleobases and amino acids; P=Psat")
# includes activities of basis species
tb <- thermo$basis   
# exclude those that are on the axes
tb <- tb[!((rownames(tb) %in% c("e-", "H2O"))),]
dp <- describe.property(c("T", "T"), c(25, 100))
db <- describe.basis(tb)
legend("bottomleft", lty=c(1, 1, NA, NA, NA), col=c("black", "red"), legend=c(dp, db))
