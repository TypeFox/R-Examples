PlotGlobEdgeCorProfile<-function(EdgeCorScore){
     N<-length(EdgeCorScore)
	 obs.s2n<-EdgeCorScore
	 obs.s2n<-sort(obs.s2n, decreasing=T)   
     location <- 1:N
     max.corr <- max(obs.s2n)
     min.corr <- min(obs.s2n)
     x <- plot(location, obs.s2n, ylab = "Correlation of edge", xlab = "Edge List Location", main = "Edge List Correlation  Profile", type = "l", lwd = 2, cex = 0.9, col = 1)    
     for (i in seq(1, N, 20)) {
       lines(c(i, i), c(0, obs.s2n[i]), lwd = 3, cex = 0.9, col = colors()[12]) # shading of correlation plot
     }
     x <- points(location, obs.s2n, type = "l", lwd = 2, cex = 0.9, col = 1)            
     lines(c(1, N), c(0, 0), lwd = 2, lty = 1, cex = 0.9, col = 1) # zero correlation horizontal line
     temp <- order(abs(obs.s2n), decreasing=T)
     arg.correl <- temp[N]
     lines(c(arg.correl, arg.correl), c(min.corr, 0.7*max.corr), lwd = 2, lty = 3, cex = 0.9, col = 1) # zero correlation vertical line

     area.bias <- signif(100*(sum(obs.s2n[1:arg.correl]) + sum(obs.s2n[arg.correl:N]))/sum(abs(obs.s2n[1:N])), digits=3)
	 phen1<-"Gain-of-correlation"
	 phen2<-"Loss-of-correlation"
     area.phen <- ifelse(area.bias >= 0, phen1, phen2)
	 delta.string <- paste("Corr. Area Bias to \"", area.phen, "\" =", 100*abs(area.bias), "%", sep="", collapse="")
     zero.crossing.string <- paste("Zero Crossing at location ", arg.correl, " (",  signif(100*arg.correl/N, digits=3), "%)")
	 leg.txt <- c(delta.string, zero.crossing.string)
     legend(x=N/10, y=max.corr, bty="n", bg = "white", legend=leg.txt, cex = 0.9)

     leg.txt <- paste("\"", phen1, "\" ", sep="", collapse="")
     text(x=1, y=-0.05*max.corr, adj = c(0, 1), labels=leg.txt, cex = 0.9)

     leg.txt <- paste("\"", phen2, "\" ", sep="", collapse="")
     text(x=N, y=0.05*max.corr, adj = c(1, 0), labels=leg.txt, cex = 0.9)
}	 