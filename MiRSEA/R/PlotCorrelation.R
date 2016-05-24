     #Results<-MsReport(MsNAME="KEGG_ERBB_SIGNALING_PATHWAY", input.ds, input.cls, weighted.score.type = 1)
	 #miRlist<-Results[[2]]
	 PlotCorrelation <-function(miRlist){ 
         miResult<-miRlist
	     N<- as.numeric(miResult[1])
		 obs.s2n<-as.numeric(miResult[3][[1]])
		 phen1<-as.character(miResult[8])
		 phen2<-as.character(miResult[9])
	 location <- 1:N
     max.corr <- max(obs.s2n)
     min.corr <- min(obs.s2n)

     x <- plot(location, obs.s2n, ylab = "differential weighted score(dw-score)", xlab = "miR List Location", main = "miR List Correlation (dw-score) Profile", type = "l", lwd = 2, cex = 0.9, col = 1)            
     for (s in seq(1, N, 20)) {
       lines(c(s, s), c(0, obs.s2n[s]), lwd = 3, cex = 0.9, col = colors()[12]) # shading of correlation plot
     }
     x <- points(location, obs.s2n, type = "l", lwd = 2, cex = 0.9, col = 1)            
     lines(c(1, N), c(0, 0), lwd = 2, lty = 1, cex = 0.9, col = 1) # zero correlation horizontal line
     temp <- order(abs(obs.s2n), decreasing=TRUE)
     arg.correl <- temp[N]
     lines(c(arg.correl, arg.correl), c(min.corr, 0.7*max.corr), lwd = 2, lty = 3, cex = 0.9, col = 1) # zero correlation vertical line

     area.bias <- signif(100*(sum(obs.s2n[1:arg.correl]) + sum(obs.s2n[arg.correl:N]))/sum(abs(obs.s2n[1:N])), digits=3)
     area.phen <- ifelse(area.bias >= 0, phen1, phen2)
     delta.string <- paste("Corr. Area Bias to \"", area.phen, "\" =", abs(area.bias), "\\%", sep="", collapse="")
     zero.crossing.string <- paste("Zero Crossing at location ", arg.correl, " (",  signif(100*arg.correl/N, digits=3), " \\%)")
     leg.txt <- c(delta.string, zero.crossing.string)
     legend(x=N/10, y=max.corr, bty="n", bg = "white", legend=leg.txt, cex = 0.9)

     leg.txt <- paste("\"", phen1, "\" ", sep="", collapse="")
     text(x=1, y=-0.05*max.corr, adj = c(0, 1), labels=leg.txt, cex = 0.9)

     leg.txt <- paste("\"", phen2, "\" ", sep="", collapse="")
     text(x=N, y=0.05*max.corr, adj = c(1, 0), labels=leg.txt, cex = 0.9)
	 }