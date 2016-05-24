file <- "/Users/Stoffi/R_files/BEAST/gigbeastcon"
file  <- "/Users/Stoffi/R_files/BEAST/cornell_strict/oxalidales_strict.tre"

tr <- read.nexus(file)


tip

tr <- read.nexus(file)



tr <- read.beast(file)
pdf("test.pdf", paper = "a4", height= 12)
plot(tr, cex = 0.7, no.margin = TRUE)
#
nodelabels(cex = 0.5)
dev.off()
system("open test.pdf")