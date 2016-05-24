## Demo for revisit(): CV of equilibrium activities of proteins in
## Pelagibacter ubique.
## Also shows use of grep.file(), read.fasta(), add.protein()
f <- system.file("extdata/fasta/HTCC1062.faa.xz", package="CHNOSZ")
# what proteins to select (set to "" for all proteins)
w <- "ribosomal"
# locate entries whose names contain w
j <- grep.file(f, w)
# get the amino acid compositions of these proteins
aa <- read.fasta(f, j)
# add these proteins to CHNOSZ's inventory
ip <- add.protein(aa)
# set up a the chemical system
basis("CHNOS+")
# calculate affinities of formation in logfO2 space
a <- affinity(O2=c(-85, -60), iprotein=ip)
# show the equilibrium activities
opar <- par(mfrow=c(2, 2))
e <- equilibrate(a, loga.balance=0, normalize=TRUE)
diagram(e, names=NULL)
# make a title
expr <- as.expression(substitute(x~y~"proteins in"~
  italic("P. ubique"), list(x=length(j), y=w)))
mtitle(c("Equilibrium activities of", expr))
# show the coefficient of variation
revisit(e, "CV")
mtitle(c("CV of equilibrium activities of", expr))
# calculate affinities in logfO2-logaH2O space
a <- affinity(O2=c(-85,-65), H2O=c(-10,0), iprotein=ip)
# show the predominances
diagram(a, normalize=TRUE, fill="heat")
# calculate the equilibrium activities
e <- equilibrate(a, loga.balance=0, normalize=TRUE)
# show the coefficient of variation
r <- revisit(e, "CV")
stopifnot(r$ix==37)
stopifnot(r$iy==53)
mtitle(c("CV of equilibrium activities of", expr))
par(opar)
