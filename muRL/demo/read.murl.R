## Specify path to .csv database of sample addresses
fpath <- system.file("data", "murljobs.csv", package = "muRL")

murljobs <- read.murl(fpath)