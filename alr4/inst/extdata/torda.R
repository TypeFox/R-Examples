# save all alr data in .rda files
library(alr3)
files <-  c(
"shocks", "turkey", "ufc", "longshoots", "caution", "sleep1", "brains",
"downer", "walleye", "MWwords", "drugcost", "galapagos", "BGSgirls", "florida",
"BGSall", "oldfaith", "lakemary", "wblake2", "heights", "mile", "water", "UN2",
"cloud", "dwaste", "physics1", "npdata", "sniffer", "wm4", "twins", "wm3",
"salarygov", "landrent", "blowAPB", "swan96", "allshoots", "wm1", "snowgeese",
"ufcwc", "jevons", "UN3", "chloride", "fuel2001", "segreg", "pipeline", "BGSboys",
"physics", "domedata", "snake", "rat", "htwt", "salary", "shortshoots", "UN1",
"challeng", "stopping", "cakes", "banknote", "highway", "wool", "Mitchell",
"galtonpeas", "lakes", "transact", "cathedral", "blowBF", "baeskel", "wm2",
"forbes", "prodscore", "turk0", "blowdown", "lathe1", "ais", "ftcollinssnow",
"hooker", "donner", "BigMac2003", "titanic", "mantel", "domedata1", "ufcdf",
"wblake", "ufcgf")

out1 <- NULL
for(f in files) {
  out1 <- c(out1,
    paste("save(", f, ", file='", f, ".rda')", sep=""))
  }
write.table(out1, "out1", row.names=FALSE, col.names=FALSE, quote=FALSE)
source("out1")
#wm5 <- read.table("http://www.stat.umn.edu/alr/data/wm5.txt", header=TRUE)
#save(wm5, file="wm5.rda")
