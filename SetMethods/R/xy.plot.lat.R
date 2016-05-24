xy.plot.lat <-
function(x, y, 
						ylim = c(-.05, 1.05), xlim = c(-.05, 1.05),
						main = "",
						pch = 19, 
						col = "black", 
						cex.fit = 1,
						ylab = "Outcome", xlab = "Condition", 
						pos.fit = "top",
						strip.cex = 0.8, 
						necessity = FALSE, 
						show.fit = TRUE, 
						case.lab = FALSE, 
						lab.pos = 4,
						labs = NULL, 
						show.hv = TRUE){

    if(necessity == TRUE){   
    # Necessity
	con <- sum(pmin(x, y))/sum(y)
	cov <- sum(pmin(x, y))/sum(x)
	cons <- format(con, digits = 3)
	storage.mode(cons) <- "numeric"
	cove <- format(cov, digits = 3)
	storage.mode(cove) <- "numeric"
    pof <- sprintf("Consistency Necessary Cond.: %.3f - Coverage Necessary Cond.: %.3f", con, cov)
	}
	else{
	# Sufficiency
	con <- sum(pmin(x, y))/sum(x)
	cov <- sum(pmin(x, y))/sum(y)
	cons <- format(con, digits = 3)
	storage.mode(cons) <- "numeric"
	cove <- format(cov, digits = 3)
	storage.mode(cove) <- "numeric"
    pof <- sprintf("Consistency Sufficient Cond.: %.3f - Coverage Sufficient Cond.: %.3f", con, cov)
	}

    
    if(show.fit == TRUE){
    	
    	if(pos.fit == "top"){
    		
    		xyplot(y ~ x | pof, 
	   	   	   ylim = ylim, 
	   	   	   xlim = xlim,
	   	   	   main = main,
	   	   	   pch = pch, 
	   	   	   col = col, 
	   	   	   ylab = ylab, 
	   	   	   xlab = xlab,
	   	   	   strip = strip.custom(par.strip.text = list(cex = strip.cex)),
	   	   	   par.settings = list(par.strip.text = list(cex = strip.cex),
	   	   	   						strip.background = list(col = NA),
	   					   	   	   layout.heights = list(strip = 1.5)),
	   					   	   
	   	   scales = list(x = list(at = seq(0, 1, .1)),
	   				 	 y = list(at = seq(0, 1, .1))), 
	   				 	 
	   		panel = function(x, y, ...){
	   			
	   			panel.abline(0, 1)
	   			panel.xyplot(x, y, ...)
	   			
	   			if(show.hv == TRUE){
	   			panel.abline(h = .5, lty = 2) 
	   			panel.abline(v = .5, lty = 2)
	   			} # end
	   			
	   			if(case.lab == TRUE){
	   			panel.text(x, y, labels = labs, pos = lab.pos)
	   			} # end
	   }# end
	   ) # end
	     
	   } # end of "top" 
	   
	   else{
    	
    	xyplot(y ~ x, 
	   	   	   ylim = ylim, xlim = xlim,
	   	   	   main = main,
	   	   	   pch = pch, col = col, 
	   	   	   ylab = ylab, xlab = xlab,
	   	   	   strip = strip.custom(par.strip.text = list(cex = strip.cex)),
	   	   	   par.settings = list(par.strip.text = list(cex = strip.cex),
	   	   	   						strip.background = list(col = NA),
	   					   	   	   layout.heights = list(strip = 1.5)),
	   					   	   
	   	   scales = list(x = list(at = seq(0, 1, .1)),
	   				 	 y = list(at = seq(0, 1, .1))), 
	   				 	 
	   		panel = function(x, y, ...){
	   			panel.abline(0, 1)
	   			panel.xyplot(x, y, ...)
	   			if(show.hv == TRUE){
	   			panel.abline(h = .5, lty = 2) 
	   			panel.abline(v = .5, lty = 2)
	   			}
	   			if(case.lab == TRUE){
	   			panel.text(x, y, labels = labs, pos = lab.pos)
	   			}
	   			panel.text(-.02, 1.02, cons, cex = cex.fit, adj = 0)
	   			panel.text(1.02, -.02, cove, cex = cex.fit, adj = 1)
	   			
	   })	   
	   } # end of corner
	   } # end of pos.fit
	   
	   else{xyplot(y ~ x, 
	   	   		   ylim = ylim, xlim = xlim,
	   	   		   main = main,
	   	   		   pch = pch, col = col, 
	   	   		   ylab = ylab, xlab = xlab,
	   	   		   strip = strip.custom(par.strip.text = list(cex = 1)),
	   	   		   par.settings = list(strip.background = list(col = NA),
	   					   	   		  layout.heights = list(strip = 1.5)),
	   			   scales = list(x = list(at = seq(0, 1, .1)),
	   				 	 		 y = list(at = seq(0, 1, .1))), 
	   				 	 
	   		panel = function(x, y, ...){
	   			panel.abline(0, 1)
	   			panel.xyplot(x, y, ...)
	   			panel.abline(h = .5, lty = 2) 
	   			panel.abline(v = .5, lty = 2)
	   			panel.text(x, y, labels = 1:40, pos = lab.pos)
	   })
	   }
}
