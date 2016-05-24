library(phylobase)
library(ape)
data(geospiza)

gtree <- extractTree(geospiza)
stopifnot(identical(gtree,prune(gtree,character(0))))

stopifnot(identical(tdata(subset(geospiza)),
                    tdata(subset(geospiza, tipLabels(geospiza)))))


tr <- read.tree(text="(((spA:0.2,(spB:0.1,spC:0.1):0.15):0.5,spD:0.7):0.2,spE:1):0.4;")
phyd <- as(tr, "phylo4d")
tipData(phyd) <- 1:5
stopifnot(identical(phyd@data,subset(phyd,tipLabels(phyd))@data))

