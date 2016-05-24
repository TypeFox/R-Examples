library(phyloclim)
library(adehabitatMA) # depends: sp
setwd("/Users/Stoffi/R_files/packages/phyloclim_archive")

## input object generation cf. "/data/samples.R"
## ---------------------------------------------
p <- samples
env <- preds
n <- 10
app <- "/Applications/maxent_3.3.3k/maxent.jar"
DIR  <- "NET"
ODIR <- paste(DIR, "out/", sep = "/")

net <- niche.equivalency.test(p = samples, env = preds, n = 9, app = app)
save(net, file = "net9rep.rda")
bst <- bg.similarity.test(p = samples, env = preds, n = 99, app = app)
save(bst, file = "bst.rda")

#load(file = "net.rda")
