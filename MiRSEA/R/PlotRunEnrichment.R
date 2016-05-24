#Results<-MsReport(MsNAME="KEGG_ERBB_SIGNALING_PATHWAY", input.ds, input.cls, weighted.score.type = 1) 
 #miRlist<-Results[[2]]
PlotRunEnrichment <-function(miRlist){ 
             miResult<-miRlist
             N<- as.numeric( miResult[1])
			 Obs.RES<-as.numeric(miResult[2][[1]])
			 obs.s2n<-as.numeric(miResult[3][[1]])
			 Obs.ES<-as.numeric(miResult[4][[1]])
			 size.M<-as.numeric(miResult[5][[1]])
			 Obs.arg.ES<-as.numeric(miResult[6][[1]])
			 Obs.indicator<-as.numeric(miResult[7][[1]])
			 phen1<-as.character(miResult[8])
		     phen2<-as.character(miResult[9])
			 MsNAME<-as.character(miResult[12])
            ind <- 1:N
            min.RES <- min(Obs.RES)
            max.RES <- max(Obs.RES)
            if (max.RES < 0.3) max.RES <- 0.3
            if (min.RES > -0.3) min.RES <- -0.3
            delta <- (max.RES - min.RES)*0.50
            min.plot <- min.RES - 2*delta
            max.plot <- max.RES
            max.corr <- max(obs.s2n)
            min.corr <- min(obs.s2n)
            Obs.correl.vector.norm <- (obs.s2n - min.corr)/(max.corr - min.corr)*1.25*delta + min.plot
            zero.corr.line <- (- min.corr/(max.corr - min.corr))*1.25*delta + min.plot
            col <- ifelse(Obs.ES > 0, 2, 4)

            sub.string <- paste("Number of miRs: ", N, " (in list), ", size.M, " (in miR set)", sep = "", collapse="")
            
            main.string <- MsNAME
            plot(ind, Obs.RES, main = main.string, sub = sub.string, xlab = "MiR List Index", ylab = "Running Enrichment Score (RES)", xlim=c(1, N), ylim=c(min.plot, max.plot), type = "l", lwd = 2, cex = 1, col = col)
            for (j in seq(1, N, 20)) {
                lines(c(j, j), c(zero.corr.line, Obs.correl.vector.norm[j]), lwd = 1, cex = 1, col = colors()[12]) # shading of correlation plot
            }
            lines(c(1, N), c(0, 0), lwd = 1, lty = 2, cex = 1, col = 1) # zero RES line
            lines(c(Obs.arg.ES, Obs.arg.ES), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = col) # max enrichment vertical line
            for (j in 1:N) {
               if (Obs.indicator[j] == 1) {
                  lines(c(j, j), c(min.plot + 1.25*delta, min.plot + 1.75*delta), lwd = 1, lty = 1, cex = 1, col = 1)  # enrichment tags
               }
            }
            lines(ind, Obs.correl.vector.norm, type = "l", lwd = 1, cex = 1, col = 1)
            lines(c(1, N), c(zero.corr.line, zero.corr.line), lwd = 1, lty = 1, cex = 1, col = 1) # zero correlation horizontal line
            temp <- order(abs(obs.s2n), decreasing=TRUE)
            arg.correl <- temp[N]
            lines(c(arg.correl, arg.correl), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = 3) # zero crossing correlation vertical line

            leg.txt <- paste("\"", phen1, "\" ", sep="", collapse="")
            text(x=1, y=min.plot, adj = c(0, 0), labels=leg.txt, cex = 1.0)

            leg.txt <- paste("\"", phen2, "\" ", sep="", collapse="")
            text(x=N, y=min.plot, adj = c(1, 0), labels=leg.txt, cex = 1.0)

            adjx <- ifelse(Obs.ES> 0, 0, 1)
           
            leg.txt <- paste("Peak at ", Obs.arg.ES, sep="", collapse="")
            text(x=Obs.arg.ES, y=min.plot + 1.8*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)

            leg.txt <- paste("Zero crossing at ", arg.correl, sep="", collapse="")
            text(x=arg.correl, y=min.plot + 1.95*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)
			}