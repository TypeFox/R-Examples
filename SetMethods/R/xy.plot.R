xy.plot <-
function(x, y, 
					ylim = c(-.05, 1.05), xlim = c(-.05, 1.05), 
					pch = 19, col = "black",
					main = "XY plot", ylab = "Outcome", xlab = "Condition",
					mar = c(4, 4, 4, 1), mgp = c(2.2, .8, 0),
					cex.fit = .6, cex.axis = .7, cex.main = 1,
					necessity = FALSE, 
					show.hv = TRUE, show.fit = TRUE, pos.fit = "top",
					case.lab = TRUE, labs = NULL, cex.lab = .8, 
					offset.x = 0, offset.y = 0, 
					pos = 4, srt = 0,
					ident = FALSE){	
						
	par(mar = mar, mgp = mgp)
	plot(x, y, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, 
		 axes = FALSE, pch = pch, main = main, cex.main = cex.main)
		 
	axis(1, at = seq(0, 1, .1), labels = seq(0, 1, .1), cex.axis = cex.axis)
	axis(2, at = seq(0, 1, .1), labels = seq(0, 1, .1), cex.axis = cex.axis, las=2)
	box(); abline(0, 1); 
	
	if(show.hv == TRUE){
		abline(v = .5, lty = 2); abline(h = .5, lty = 2)
		}
	
   
    if(necessity == TRUE){   
    # Necessity
	con <- sum(pmin(x, y))/sum(y)
	cov <- sum(pmin(x, y))/sum(x)
	cons <- format(con, digits = 3)
	storage.mode(cons) <- "numeric"
	cove <- format(cov, digits = 3)
	storage.mode(cove) <- "numeric"
    lab <- sprintf("Consistency Necessary Condition: %.3f - Coverage Necessary Condition: %.3f", 
    con, cov)
    cons.c <- paste("Consistency Necessary Condition", cons, sep = ": ")
    cove.c <- paste("Coverage Necessary Condition", cove, sep = ": ")

	}
	else{
	# Sufficiency
	con <- sum(pmin(x, y))/sum(x)
	cov <- sum(pmin(x, y))/sum(y)
	cons <- format(con, digits = 3)
	storage.mode(cons) <- "numeric"
	cove <- format(cov, digits = 3)
	storage.mode(cove) <- "numeric"
    lab <- sprintf("Consistency Sufficient Condition: %.3f - Coverage Sufficient Condition: %.3f", 
    con, cov)
	}
	
    if(show.fit == TRUE){
    	
    	if(pos.fit == "top" ){	
    	mtext(lab, line = 0.3, cex = cex.fit)
    	}
    	
    	if(pos.fit == "corner" ){
    		text(-.05, 1.05, cons, cex = cex.fit, adj = 0)
    		text(1.05, -.05, cove, cex = cex.fit, adj = 1)
    	}
		}
    	
    if(case.lab == TRUE){
    	text(x + offset.x, y + offset.y, 
    		 labels = labs, cex = cex.lab, pos = pos, srt = srt)
    	} 
    	     
    if(ident == TRUE){
    	id <- identify(x, y, labels = labs, cex = cex.lab)
    	}
}
