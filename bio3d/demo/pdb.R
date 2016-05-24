###
### Example of PDB file manipulation, searching, alignment etc.
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
## Basic PDB file reading and manipulation  #
##                                          #
#############################################
pause()

# Read an online RCSB Protein Data Bank structure
pdb <- read.pdb("4q21")

# Whats in the new pdb object
print(pdb)

pause()

# Most bio3d functions, including read.pdb(), return list objects
attributes(pdb)

pdb$atom[1:3, c("resno", "resid", "elety", "x", "y", "z")]

pause()

# Selection of substructure regions with 'atom.select()'' function
inds <- atom.select(pdb, elety = c("N","CA","C"), resno=4:6)

pdb$atom[inds$atom,]

pause()

# Simple B-factor plot
ca.inds <- atom.select(pdb, "calpha")
plot.bio3d( pdb$atom[ca.inds$atom,"b"], sse=pdb, ylab="B-factor")


###################################
##                                #
## Search for similar structures  #
##                                #
###################################

# Use sequence
aa <- pdbseq(pdb)
aa

pause()

# Blast the RCSB PDB to find similar sequences 
blast <- blast.pdb(aa)
head(blast$hit.tbl)

pause()

# Plot results
top.hits <- plot(blast)
head(top.hits$hits)

pause()

## Download and and analyze further ....
## raw.files <- get.pdb(top.hits$pdb.id, path="raw_hits")
## files <- pdbsplit(raw.files, top.hits$hits, path="top_hits")
## pdbs <- pdbaln(files)
## ...<ETC>...
