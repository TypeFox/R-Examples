require(lattice)
require(SoDA)
## requires mars data.frame
SoDA:::.getMarsData()
pdf(file = "Examples/marsPlot.pdf", width = 7.5, height = 3)
par(cex=.6)
plot(mars$Date, mars$Declination)
dev.off()

