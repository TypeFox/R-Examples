###
### Example of PCA on a collection of PKA structures
###    and a large collection of transducin structure
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

################################################
##                                             #
## Basic PCA of related X-ray structures       #
## (requires the 'muscle' program installed)   #
##                                             #
################################################

pause()

### Set temp dir to store PDB files
tmp.dir <- tempdir()

## Specify PDB identifiers
ids <- c("1cdk_A", "3agm_A", "1cmk_E",
         "3dnd_A", "1q8w_A")

## Download PDBs
raw.files <- get.pdb(ids, path=tmp.dir)

pause()

## Split PDBs by chain ID
files <- pdbsplit(raw.files, ids, path=tmp.dir)

pause()

## Sequence/structure alignment
pdbs <- pdbaln(files)

pause()

## Find invariant core
core <- core.find(pdbs)

pause()

## Fit structures to core region
xyz <- pdbfit(pdbs, inds=core$c1A.xyz)
           ## outpath="core_fit/", full.pdbs=T, het2atom=T)

pause()

## Locate gap containing positions
gaps.pos <- gap.inspect(pdbs$xyz)

## Perform PCA on non-gap containing positions
pc.xray <- pca.xyz(xyz[,gaps.pos$f.inds])

pause()

## Plot x-ray results
plot(pc.xray)

pause()

#############################################
##                                          #
## Larger transducin example                #
##                                          #
#############################################

data(transducin)
attach(transducin, warn.conflicts=FALSE)

## data 'transducin' contains objects
## - pdbs:       aligned C-alpha coordinates for 53 transducin 
##               structures from the PDB
## - annotation: annotation of the 53 PDBs

## Note that this data can be generated from scratch by following the
## Comparative Structure Analysis with Bio3D Vignette available both
## on-line and from within the Bio3D package.

pdbs <- transducin$pdbs
annotation <- transducin$annotation

pause()

## Inspect gaps
gaps.pos <- gap.inspect(pdbs$xyz)

## Previously fitted coordinates invariance core
xyz <- pdbs$xyz

## Do PCA
pc.xray <- pca.xyz(xyz[, gaps.pos$f.inds])

pause()


## Plot overview
plot(pc.xray, col=annotation[, "color"])

## Plot atom wise loadings
plot.bio3d(pc.xray$au[,1], ylab="PC1 (A)")

pause()

unlink(tmp.dir)






