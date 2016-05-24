## From ~/soar/attie_mouse
## Phase 0.
## Run-specific settings:
## Must specify "tissue" and "runnum" (and "subset.sex" if used) before sourcing this file.

library(B6BTBR07n)
library(qtlhot)

offset <- 0
if(!exists("subset.sex")) subset.sex <- NULL

Nmax <- 2000
drop.lod <- 1.5
size.set <- 250
datadir <- "~/soar/diabetes48" ## Aimee's diabetes4 folder.

cat("link needed data\n")
## tissue.dir/control_probes.txt
tissue.dir <- ifelse(datadir == ".", ".", file.path(datadir, "final.Rdata"))

## lod.thr is in null.dir/quantiles.csv
null.dir <- ifelse(datadir == ".", ".", file.path(datadir, "nullsims"))
system(paste("ln -s", file.path(null.dir, "quantiles.csv"), "."))
lod.thrs <- read.csv(file.path(null.dir, "quantiles.csv"), header = TRUE)
tmp <- substring(lod.thrs[,1], 1, 2)
if(is.null(subset.sex)) {
  lod.thrs <- rev(sort(lod.thrs$full))
} else {
  lod.thrs <- rev(sort(lod.thrs[[tolower(subset.sex)]]))
  if(tolower(subset.sex) == "full")
    subset.sex <- NULL
}
names(lod.thrs) <- rev(tmp)
  
## mlratio.file = tissue.dir/tissue_mlratio_final.RData
## contains tissue.mlratio data.frame.
trait.file <- paste(tissue, "mlratio", "final.RData", sep = "_")
trait.mlratio <- file.path(tissue.dir, trait.file)
trait.file <- "trait.RData"
system(paste("ln -s", trait.mlratio, trait.file))
trait.matrix <- paste(tissue, "mlratio", sep = ".")
trait.index <- "MouseNum"

## Control probes. Not needed any more.
## droptrait.file <- "control_probes.txt"
## system(paste("ln -s", file.path(tissue.dir, droptrait.file), "."))
## droptrait.names <- read.table(file.path(tissue.dir, droptrait.file), header = FALSE)[,1]
droptrait.names <- NULL
## Get cross with Sex and batch.
sex <- "Sex"
trait.index <- "MouseNum"
batch.effect <- paste(tissue, "batch", sep = ".")
  
########################################################################################
## Annotation
cat("annotation\n")
keeptrait.names <- B6BTBR07n.annotation$a_gene_id[B6BTBR07n.annotation$probe_use >= 0]
cat("big phase0\n")
big.phase0(".", B6BTBR07n, trait.file, trait.matrix, droptrait.names,
           keeptrait.names,
           lod.thrs,
           sex = sex, trait.index = trait.index,
           batch.effect = batch.effect, size.set = size.set,
           offset = offset, subset.sex = subset.sex, verbose = TRUE)

########################################################################################
# Need to set runnum, tissue.
cat("make local SOAR directories\n")
run <- paste(tissue, runnum, sep = "")

system(paste("mkdir", run))
system(paste("mkdir ", run, "/shared", sep = ""))
system("rm -f Trait.RData")
system(paste("mv *.RData ", run, "/shared", sep = ""))
system(paste("sed -e 's/TISSUE/", tissue, "/' < ../params.txt > ", run, "/params.txt", sep = ""))

cat("submit SOAR job to simon.stat\n")
## must be on simon.stat for these:
system(paste("cp -r ", run, "~/soar/rundata/stage_qtlhot"))
system(paste("cd ~/soar/rundata/stage_qtlhot; ./stage.pl --input", run))
