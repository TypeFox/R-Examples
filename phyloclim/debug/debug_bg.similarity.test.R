library(phyloclim)
setwd("/Users/Stoffi/R_files/packages/phyloclim_archive")

maxent.exe <- paste(system.file(package="dismo"), 
                    "/java/maxent.jar", sep = "")
#spec <- c("adenophylla", "arenaria")
species <- c("enneaphylla", "laciniata")
data(sites)
samples <- sites[grep(paste(species, collapse = "|"), sites$spec), ]
data.path <- system.file("extdata", package = "phyloclim")
preds <- list.files(path = data.path, pattern = "[.]asc")
preds <- paste(data.path, preds, sep = "/")
preds <- stack(lapply(X = preds, FUN = raster))
reps <- n <- 9
env <- preds; p <- samples
bst.ennlac.95 <- bg.similarity.test(samples, preds, reps, .95, maxent.exe, dir = "bst95")
bst.ennlac.99 <- bg.similarity.test(samples, preds, reps, .99, maxent.exe)
net.ennlac <- niche.equivalency.test(samples, preds, reps, maxent.exe, "net")
