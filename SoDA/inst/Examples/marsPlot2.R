## requires mars data.frame
SoDA:::.getMarsData()
pdf(file = "Examples/marsPlot2.pdf", width = 7.5, height = 3)
par(cex=.6)
noTime <- nchar(mars$Time)==0
pchTime <- ifelse(noTime, 3, 1)
plot(mars$Date, mars$Declination,
      xlab = "Date of Observation",
     ylab =  "Declination of Mars",
     pch = pchTime)
# mtext("+", at = Date[noTime], cex = .5)
dev.off()

