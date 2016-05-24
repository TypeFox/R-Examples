plot.HTSClusterWrapper <-
function (x, file.name = FALSE, 
graphs = c("capushe", "ICL", "BIC"), capushe.validation=NA, ...) 
{	
    	if (class(x) != "HTSClusterWrapper") {
        	stop(paste(sQuote("x"), sep = ""), " must be of class ", 
            paste(dQuote("HTSClusterWrapper"), sep = ""), sep = "")
    	}
        
	if(file.name != FALSE) pdf(paste(file.name));
        
        if("ICL" %in% graphs & "BIC" %in% graphs) {
          par(mfrow = c(1,2), mar = c(4,4,2,2))
          gpl <- unlist(lapply(strsplit(names(x$ICL.all), split="=", fixed=TRUE), function(xx) xx[2]))
          plot(gpl, x$ICL.all, xlab = "Number of clusters", ylab = "ICL", main="ICL", pch=19)
          lines(gpl, x$ICL.all, lwd=2)
          points(gpl[which(x$ICL.all == max(x$ICL.all))], max(x$ICL.all), col="red", pch="X", font=2, cex=2)
          plot(gpl, unlist(lapply(x$all.results, function(xx) xx$BIC)), xlab = "Number of clusters", ylab = "BIC",
               main="BIC", pch=19)
          lines(gpl, unlist(lapply(x$all.results, function(xx) xx$BIC)), lwd=2)          
        }

        if("ICL" %in% graphs & !"BIC" %in% graphs) {
          par(mar = c(4,4,2,2))
          gpl <- unlist(lapply(strsplit(names(x$ICL.all), split="=", fixed=TRUE), function(xx) xx[2]))
          plot(gpl, x$ICL.all, xlab = "Number of clusters", ylab = "ICL", main="ICL", pch=19)
          lines(gpl, x$ICL.all, lwd=2)
          points(gpl[which(x$ICL.all == max(x$ICL.all))], max(x$ICL.all), col="red", pch="X", font=2, cex=2)
        }

        if(!"ICL" %in% graphs & "BIC" %in% graphs) {
          par(mar = c(4,4,2,2))
          gpl <- unlist(lapply(strsplit(names(x$ICL.all), split="=", fixed=TRUE), function(xx) xx[2]))
          plot(gpl, unlist(lapply(x$all.results, function(xx) xx$BIC)), xlab = "Number of clusters", ylab = "BIC",
               main="BIC", pch=19)
          lines(gpl, unlist(lapply(x$all.results, function(xx) xx$BIC)), lwd=2)          
        }
		
		if("capushe" %in% graphs) {
			plot(x$capushe, newwindow=FALSE)
			if(is.na(capushe.validation) == FALSE) {
				Kchoice <- as.numeric(unlist(lapply(strsplit(names(x$logLike), "="), function(x) x[2])))
				if(capushe.validation >= max(Kchoice)) stop("Number of clusters for capushe validation should be less than largest number of clusters.");
				np <- (Kchoice-1) + (length(unique(x$all.result[[1]]$conds))-1)*(Kchoice)
				n <- nrow(x$all.result[[1]]$probaPost)
				mat <- cbind(Kchoice, np/n, np/n, -x$logLike.all)
				ResCapushe <- capushe(mat, n)
				validation(ResCapushe, mat[-c(which(Kchoice < capushe.validation)),], 
					newwindow=FALSE)
			}
		}
		

	if(file.name != FALSE)  dev.off();
      }

