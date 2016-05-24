## run.wjd with proteins: cell periphery of yeast
# get the proteins in the requested location
y <- yeastgfp("cell.periphery")
# get the amino acid compositions of the proteins
aa <- more.aa(y$protein, "Sce")
# don't use those with NA abundance or sequence
ina <- is.na(y$abundance) | is.na(aa$chains)
aa <- aa[!ina, ]
# let's try normalizing the proteins to single residues
# columns 6:25 are the actual amino acid counts
aa.625 <- aa[, 6:25]
aa[, 6:25] <- aa.625 / rowSums(aa.625)
# add proteins to thermo$protein
add.protein(aa)
# add proteins to thermo$obigt
iobigt <- info(paste(aa$protein, aa$organism, sep="_"))
# use equal initial abundances, with total equal to yeastGFP abundances
Y <- rep(mean(y$abundance[!ina]), length(y$abundance[!ina]))
# run the Gibbs energy minimization
w <- run.wjd(iobigt, Y=Y, imax=100)
# make a log-log plot
plot(log10(y$abundance[!ina]), log10(w$X), xlim=c(1.5, 5), ylim=c(1.5, 5),
  xlab="log10(abundance) reported in YeastGFP study",
  ylab="log10(abundance) calculated using Gibbs energy minimization")
# get the element potentials (tolerating "close enough" to equilibrium)
emu <- equil.potentials(w, tol=1e7)
# then the logarithms of activities of the basis species
basis("CHNOS")
bl <- basis.logact(emu)
# make a title and legend
title(main="relative abundances of proteins: yeast cell periphery")
basis(names(bl), bl)
legend("topleft", describe.basis(digits=2))
