## requires mars data.frame
SoDA:::.getMarsData()
pdf(file = "Examples/marsPlot3.pdf", width = 7.5, height = 3)
xx = par("mar")
xx[[1]] = xx[[1]]-2
par(mar = xx)
par(cex=.6)
noTime <- nchar(mars$Time)==0
pchTime <- ifelse(noTime, 3, 1)
plot(mars$Date, mars$Declination,
      xlab = "Date of Observation",
     ylab =  "Declination of Mars",
     pch = pchTime)
legend("right", c("Time of Day Given", "No Time of Day"),
                        pch=c(1,3))
points(mars$Date[noTime], rep(31, sum(noTime)), xpd=TRUE, pch=3)
# mtext("+", at = Date[noTime])
dev.off()

