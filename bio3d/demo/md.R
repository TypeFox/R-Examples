###
### Example of basic molecular dynamics trajectory analysis
###
### Authors Xin-Qiu Yao
###         Lars Skjaerven
###         Barry J Grant
###
require(bio3d); require(graphics);

pause <- function() {
  cat("Press ENTER/RETURN/NEWLINE to continue.")
  readLines(n=1)
  invisible()
}

#############################################
##                                          #
## Basic analysis of HIVpr trajectory data  #
##                                          #
#############################################
pause()

# Read example trajectory file
trtfile <- system.file("examples/hivp.dcd", package="bio3d")
trj <- read.dcd(trtfile)

# Read the starting PDB file to determine atom correspondence
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
pdb <- read.pdb(pdbfile)

# Whats in the new pdb object
print(pdb)

pause()

# How many rows (frames) and columns (coords) present in trj
dim(trj)
ncol(trj) == length(pdb$xyz)
pause()

# Trajectory Frame Superposition on Calpha atoms
ca.inds <- atom.select(pdb, elety = "CA")
xyz <- fit.xyz(fixed = pdb$xyz, mobile = trj, 
	           fixed.inds = ca.inds$xyz, 
	           mobile.inds = ca.inds$xyz)

# Root Mean Square Deviation (RMSD)
rd <- rmsd(xyz[1, ca.inds$xyz], xyz[, ca.inds$xyz])
plot(rd, typ = "l", ylab = "RMSD", xlab = "Frame No.")
points(lowess(rd), typ = "l", col = "red", lty = 2, lwd = 2)

summary(rd)
pause()

# Root Mean Squared Fluctuations (RMSF)
rf <- rmsf(xyz[, ca.inds$xyz])
plot(rf, ylab = "RMSF", xlab = "Residue Position", typ="l")

pause()

# Principal Component Analysis
pc <- pca.xyz(xyz[, ca.inds$xyz])
plot(pc, col = bwr.colors(nrow(xyz)))

pause()

# Cluster in PC space
hc <- hclust(dist(pc$z[, 1:2]))
grps <- cutree(hc, k = 2)
plot(pc, col = grps)

pause()

# Cross-Correlation Analysis
cij <- dccm(xyz[, ca.inds$xyz])
plot(cij)
## view.dccm(cij, pdb, launch = TRUE)
