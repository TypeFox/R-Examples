## Oxygen fugacity - activity of H2O predominance 
## diagrams for proteologs for 23 YeastGFP localizations
# use old properties of [Met] to reproduce this example
data(thermo)
add.obigt()
# arranged by decreasing metastability:
# order of this list of locations is based on the 
# (dis)appearance of species on the current set of diagrams
names <- c("vacuole", "early.Golgi", "ER", "lipid.particle",
  "cell.periphery", "ambiguous", "Golgi", "mitochondrion",
  "bud", "actin", "cytoplasm", "late.Golgi",
  "endosome", "nucleus", "vacuolar.membrane", "punctate.composite",
  "peroxisome", "ER.to.Golgi", "nucleolus", "spindle.pole",
  "nuclear.periphery", "bud.neck", "microtubule")
nloc <- c(4, 5, 3, 4, 4, 3)
# define the system
basis("CHNOS+")
# get protein names and abundances in each location
gfp <- yeastgfp(names)
# get amino acid compositions of proteins
aa <- more.aa(gfp$protein, "Sce")
# calculate average amino acid compositions 
for(i in 1:length(names)) {
  avgaa <- aasum(aa[[i]], gfp$abundance[[i]], average=TRUE, protein=names[i])
  add.protein(avgaa)
}
species(names, "Sce")
a <- affinity(H2O=c(-5, 0, 256), O2=c(-80, -66, 256))
# setup the plot
layout(matrix(c(1, 1,2:7), byrow=TRUE, nrow=4), heights=c(0.7, 3, 3, 3))
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, paste("Subcellular proteins of S. cerevisiae,",
 "after Dick, 2009\n"), cex=1.5)
text(0.5, 0.2, describe.basis(ibasis=c(1, 3, 4, 6), oneline=TRUE), cex=1.5)
opar <- par(mar=c(3, 4, 1, 1), xpd=TRUE)
fill <- heat.colors(length(names))
inames <- 1:length(names)
for(i in 1:length(nloc)) {
  diagram(a, normalize=TRUE, names=names[inames], groups=as.list(inames),
    fill=fill[inames], cex.axis=0.75, cex.names=1)
  label.plot(letters[i], xfrac=0.95, yfrac=0.9, paren=TRUE, italic=TRUE)
  title(main=paste(length(inames), "locations"))
  # take out the stable species
  inames <- inames[-(1:nloc[i])]
}
# return to plot defaults
layout(matrix(1))
par(opar)
# reset thermodynamic database
data(thermo)
