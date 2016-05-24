library(phyloclim)
library(raster)
setwd("/Users/Stoffi/R_files/packages/phyloclim_archive")

# ENMs
# ---------------
spec <- c("loricata", "enneaphylla", "laciniata")
fn <- paste("../../maxent/out_div_cas_bio19fg/", spec,
            "_layers_div1M.asc", sep = "")
enm <- lapply(X = fn, FUN = read.asciigrid)

species <- c("enneaphylla", "laciniata")
samples <- read.csv("/Users/Stoffi/R_files/oxalis/diversity/data/oxalis_bio19f.csv")
samples <- samples[samples$spec %in% species, 1:3]

# predictor variables: mean temperature (bio1) + mean precipitation (bio11)
# -------------------------------------------------------------------------
meanT <- read.asciigrid("/Users/Stoffi/R_files/maxent/layers_div1M/bio1.asc")
meanP <- read.asciigrid("/Users/Stoffi/R_files/maxent/layers_div1M/bio12.asc")
preds <- cbind(meanT, meanP)
  
#preds <- stack(list(raster(meanT), raster(meanP)))
#layerNames(preds) <- c("mean_temperature", "mean_precipitation")
# plot(preds)

p <- samples
x <- preds
n <- 10
app <- "/Applications/maxent_3.3.3k/maxent.jar"
DIR  <- "NET"
ODIR <- paste(DIR, "out/", sep = "/")

net <- niche.equivalency.test(p = samples, x = preds, n = 15, app = app, dir = "NET")
bst <- bg.similarity.test(p = samples, x = preds, n = 40, app = app, dir = "BST")

