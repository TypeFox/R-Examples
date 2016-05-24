PlotRunEnrichment<-function(EdgeCorScore,PathwayResult,weighted.score.type = 1){
           N<-length(EdgeCorScore)
           obs.s2n<-EdgeCorScore
		   Result<-as.matrix(PathwayResult[[1]])
		   Obs.arg.ES<-as.numeric(Result[,"List Loc"])
		   Obs.ES<-as.numeric(Result[,"RES"])
		   gene.list<-order(obs.s2n, decreasing=T)
		   gene.set<-match(Result[,"EdgeID"], names(obs.s2n)) 
		   obs.s2n   <- sort(obs.s2n, decreasing=T) 
		   correl.vector<-obs.s2n
		   tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
           no.tag.indicator <- 1 - tag.indicator 
           N <- length(gene.list) 
           Nh <- length(gene.set) 
           Nm <-  N - Nh 
           if (weighted.score.type == 0) {
           correl.vector <- rep(1, N)
           }
           alpha <- weighted.score.type
           correl.vector <- abs(correl.vector**alpha)
           sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
           norm.tag    <- 1.0/sum.correl.tag
		   norm.no.tag <- 1.0/Nm
           Obs.RES<- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)   
		   		   
		   location<-length(which(Result[,"CORE_ENRICHMENT"]=="YES"))
		   Obs.arg.ES<-Obs.arg.ES[location]
		   Obs.ES<-Obs.ES[location]
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

            sub.string <- paste("Number of edges: ", N, " (in list), ", length(Result[,1]), " (in pathway)", sep = "", collapse="")
            main.string <- names(PathwayResult)
            plot(ind, Obs.RES, main = main.string, sub = sub.string, xlab = "Edge List Index", ylab = "Running Enrichment Score (RES)", xlim=c(1, N), ylim=c(min.plot, max.plot), type = "l", lwd = 2, cex = 1, col = col)
            for (j in seq(1, N, 20)) {
               lines(c(j, j), c(zero.corr.line, Obs.correl.vector.norm[j]), lwd = 1, cex = 1, col = colors()[12]) # shading of correlation plot
            }
            lines(c(1, N), c(0, 0), lwd = 1, lty = 2, cex = 1, col = 1) # zero RES line
			
            lines(c(Obs.arg.ES, Obs.arg.ES), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = col) # max enrichment vertical line
            for (j in 1:N) {
               if (is.element(j,Result[,"List Loc"])) {
                  lines(c(j, j), c(min.plot + 1.25*delta, min.plot + 1.75*delta), lwd = 0.1, lty = 1, cex = 1, col = 1)  # enrichment tags
               }
            }
            lines(ind, Obs.correl.vector.norm, type = "l", lwd = 1, cex = 1, col = 1)
            lines(c(1, N), c(zero.corr.line, zero.corr.line), lwd = 1, lty = 1, cex = 1, col = 1) # zero correlation horizontal line
            temp <- order(abs(obs.s2n), decreasing=T)
            arg.correl <- temp[N]
            lines(c(arg.correl, arg.correl), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = 3) # zero crossing correlation vertical line

            leg.txt <- paste("\"", "Gain-of-correlation", "\" ", sep="", collapse="")
            text(x=1, y=min.plot, adj = c(0, 0), labels=leg.txt, cex = 1.0)

            leg.txt <- paste("\"", "Loss-of-correlation", "\" ", sep="", collapse="")
            text(x=N, y=min.plot, adj = c(1, 0), labels=leg.txt, cex = 1.0)

            adjx <- ifelse(Obs.ES > 0, 0, 1)
           
            leg.txt <- paste("Peak at ", Obs.arg.ES, sep="", collapse="")
            text(x=Obs.arg.ES, y=min.plot + 1.8*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)

            leg.txt <- paste("Zero crossing at ", arg.correl, sep="", collapse="")
            text(x=arg.correl, y=min.plot + 1.95*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)
}