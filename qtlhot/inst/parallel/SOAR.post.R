## From ~/soar/attie_mouse
## Must be run on simon; 
## Need to set up job, tissue, runnum, phase3.done.
## job <- 203
## tissue <- "islet"
## runnum <- 1
## phase3.done <- TRUE

main <- tissue

## Kludge to allow for no phase3 vs. redo phase3
if(!exists("phase3.done")) {phase3.done <- TRUE}
if(!exists("phase3.none")) {phase3.none <- TRUE} else {phase3.done <- FALSE}

## Kludge for sex.
if(!exists("subset.sex")) subset.sex <- NULL
if(!is.null(subset.sex)) {if(tolower(subset.sex) == "full") subset.sex <- NULL}

library(qtlhot)

## Specific locations of stuff.
run <- paste(tissue, runnum, sep = "")
cross.index <- 1
datadir <- "~/soar/diabetes48"

if(phase3.done) {
  ## Can be done on any machine.
  cat("phase3 done\n")
  dirpath <- file.path("soar", "results", job, paste("dataset_", run, "-1", sep = ""))
  cat("scp from simon", dirpath, "\n")
  system(paste("scp yandell@simon:", file.path(dirpath, "Phase31.RData"), " ", run, sep = ""))
  load(file.path(run, "Phase31.RData"))
  system(paste("rm -f", file.path(run, "Phase31.RData")))
} else {
  ## Must be on simon to do this!
  cat("Running Phase 3 (may take some time) ...\n")
  if(phase3.none) {
    dirpath <- file.path("~/soar", "condorruns", "qtlhot", job, paste("dataset_", run, "-1", sep = ""))
    dirpath2 <- dirpath
  } else {
    dirpath <- file.path("~/soar", "results", "qtlhot", job, paste("dataset_", run, "-1", sep = ""))
    dirpath2 <- file.path(dirpath, "phase2data")
  }
  qtlhot.phase3(dirpath2, cross.index, dirpath.save = ".", big = TRUE, verbose = TRUE)
  load("Phase31.RData")
}

## Single trait thresholds.
## This stuff should be in Phase3, but is not yet.
null.dir <- file.path(datadir, "nullsims")
lod.thrs <- read.csv(file.path(null.dir, "quantiles.csv"), header = TRUE)
## These are actual probabilities used for quantiles, but not in that file.
probs <- seq(.8,.99, by = .01)

tmp <- substring(lod.thrs[,1], 1, 2)
if(is.null(subset.sex)) {
    lod.thrs <- rev(sort(lod.thrs$full))
  } else {
      lod.thrs <- rev(sort(lod.thrs[[tolower(subset.sex)]]))
    }
names(lod.thrs) <- rev(tmp)

###########################################################################################################
## Thresholds for different significance levels.

pdf(paste("quant", run, job, "pdf", sep = "."), width = 7, height = 7)
out <- quant.plot(max.lod.quant, max.N, max.N.window, lod.thrs, probs, level = 0.95)
dev.off()

## Why are there 2 threshold values?
sink(paste("quant", tissue, job, "txt", sep = "."))
options(width=200)
round(out$quant.N.window, 3)
round(out$quant.N, 3)
sink()

#quant.level <- out$quant.level
#quant.thr <- out$quant.thr
save(lod.thrs, out, file = paste("thr", run, job, "RData", sep = "."), compress = TRUE)

###########################################################################################################
## Landscape of hotspots with 5% threshold.

## Load unpermuted data (first set).
if(phase3.done) {
  system(paste("scp yandell@simon:",
               file.path(dirpath, "phase2data", paste("perm.", cross.index, "_1.RData", sep = "")),
               " ", run, sep = ""))
} else {
  system(paste("scp yandell@simon:",
               file.path(dirpath2, paste("perm.", cross.index, "_1.RData", sep = "")),
               " ", run, sep = ""))
}
load(file.path(run, paste("perm.", cross.index, "_1.RData", sep = "")))
load(paste("thr", run, job, "RData", sep = "."))

library(B6BTBR07n)

## Need Phase1.RData for n.traits and cross
system(paste("scp yandell@simon:", file.path(dirpath, "phase2data", "Phase1.RData"),
             " ", run, sep = ""))
attach(file.path(run, "Phase1.RData"), warn.conflicts = FALSE)
n.traits <- length(all.traits)
cross <- cross
detach()
system(paste("rm -f", file.path(run, "Phase1.RData")))

#load(file.path(run, "shared", "cross.RData"))
## This step takes awhile. And it may still be broken for subset.sex=full
cat("hotspot scan ...\n")
hotspots <- hotspot.scan(cross, scan.hl, out$lod.thr, out$quant.level, window = 5)
save(hotspots, file = paste("hot", run, job, "RData", sep = "."), compress = TRUE)

pdf(paste("hot", run, job, "pdf", sep = "."), width = 7, height = 7)
hotspot.plot(hotspots, quant.thr = out$quant.thr, maps = B6BTBR07n.maps, main = main)
dev.off()
