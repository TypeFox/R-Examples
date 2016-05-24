###
### Examples from NMA Vignette
###
### Authors Lars Skjaerven
###         Xin-Qiu Yao
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
## Basic usage                              #
##                                          #
#############################################

### Read PDB and Calculate Normal Modes
pdb <- read.pdb("1hel")
modes <- nma(pdb)

pause()

### Print a summary
print(modes)

pause()

### Plot the nma object for a quick overview
plot(modes)

pause()

### Calculate cross-correlations
cm <- dccm(modes)

pause()

### Plot correlation map
plot(cm, sse=pdb)

pause()

### Calculate modes with force field ANM
modes.anm <- nma(pdb, ff="anm")

pause()

### Investigate modes similarity with RMSIP
r <- rmsip(modes, modes.anm)

pause()

### Plot RMSIP results
plot(r, xlab="ANM", ylab="C-alpha FF")

pause()

################################################
##                                             #
## Ensemble NMA                                #
## (requires the 'muscle' program installed)   #
##                                             #
################################################

pause()

### Set temp dir to store PDB files
tmp.dir <- tempdir()

### Download a set of DHFR structures 
ids <- c("1rx2_A", "1rx4_A", "1rg7_A", 
         "3fyv_X", "3sgy_B")



### Download and split by chain ID
raw.files <- get.pdb(ids, path=tmp.dir)

pause()

### Split PDB files by chain ID
files <- pdbsplit( raw.files, ids, path=tmp.dir)

pause()

### Align structures
pdbs <- pdbaln(files)

pause()

### View sequence identity
summary( c(seqidentity(pdbs)) )

pause()

### Calculate modes of aligned proteins
modes <- nma(pdbs)

pause()

## Print a summary
print(modes)

pause()

### Plot fluctuations
plot(modes, pdbs)

pause()

### Cluster Modes simiarlity
heatmap(1-modes$rmsip, labCol=ids)

unlink(tmp.dir)
