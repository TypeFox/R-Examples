PSRplot <- function(x, info = c("area", "null", "Brownian", "both"), ...){

	if(info == "both"){
		
		par(mfrow = c(2, 1))
		
		plot(x@PSR, 
				xlab = "cumulative eigenvalues(%)", 
				ylab = "r squared", ylim = c(0,1), 
				xlim = c(0,1), bty = "L", expected = segments(0, 0, 1, 1), pch = 16, ...)	
		e <- rowMeans(x@nullPSR)
		f <- apply(x@nullPSR, 1, sd)
		f <- f/sqrt(ncol(x@nullPSR))
		coordsShade <- data.frame(x = c(x@PSR$Cumul.eigen.values, sort(x@PSR$Cumul.eigen.values, decreasing = T)), 
				y = c((e + (1.96 * f)), sort((e - (1.96 * f)), decreasing = T)))
		
		polygon(coordsShade, border = NA,  col = "#FFFF0070");
		points(x@PSR$r.squared ~ x@PSR$Cumul.eigen.values, pch = 16);
		lines(e ~ x@PSR$Cumul.eigen.values, lty = 2, cex = 1.5)
		
		plot(x@PSR, 
				xlab = "cumulative eigenvalues(%)", 
				ylab = "r squared", ylim = c(0,1), 
				xlim = c(0,1), bty = "L", expected = segments(0, 0, 1, 1), pch = 16, ...)	
		e <- rowMeans(x@BrownianPSR)
		f <- apply(x@BrownianPSR, 1, var)
		coordsShade <- data.frame(x = c(x@PSR$Cumul.eigen.values, sort(x@PSR$Cumul.eigen.values, decreasing = T)), 
				y = c((e+f), sort((e-f), decreasing = T)))
		
		polygon(coordsShade, border = NA,  col = "#FF000070");
		points(x@PSR$r.squared ~ x@PSR$Cumul.eigen.values, pch = 16);
		lines(e ~ x@PSR$Cumul.eigen.values, lty = 2, cex = 1.5)
		
		par(mfrow = c(1, 1))
	} else {
		
		if(info == "null"){
			
			plot(x@PSR, 
					xlab = "cumulative eigenvalues(%)", 
					ylab = "r squared", ylim = c(0,1), 
					xlim = c(0,1), bty = "L", expected = segments(0, 0, 1, 1), pch = 16, ...)	
			e <- rowMeans(x@nullPSR)
			f <- apply(x@nullPSR, 1, sd)
			f <- f/sqrt(ncol(x@nullPSR))
			coordsShade <- data.frame(x = c(x@PSR$Cumul.eigen.values, sort(x@PSR$Cumul.eigen.values, decreasing = T)), 
					y = c((e + (1.96 * f)), sort((e - (1.96 * f)), decreasing = T)))
			
			polygon(coordsShade, border = NA,  col = "#FFFF0070");
			points(x@PSR$r.squared ~ x@PSR$Cumul.eigen.values, pch = 16);
			lines(e ~ x@PSR$Cumul.eigen.values, lty = 2, cex = 1.5)
		}	
		
		if(info == "Brownian"){
			
			plot(x@PSR, 
					xlab = "cumulative eigenvalues(%)", 
					ylab = "r squared", ylim = c(0,1), 
					xlim = c(0,1), bty = "L", expected = segments(0, 0, 1, 1), pch = 16, ...)	
			e <- rowMeans(x@BrownianPSR)
			f <- apply(x@BrownianPSR, 1, var)
			coordsShade <- data.frame(x = c(x@PSR$Cumul.eigen.values, sort(x@PSR$Cumul.eigen.values, decreasing = T)), 
					y = c((e+f), sort((e-f), decreasing = T)))
			
			polygon(coordsShade, border = NA,  col = "#FF000070");
			points(x@PSR$r.squared ~ x@PSR$Cumul.eigen.values, pch = 16);
			lines(e ~ x@PSR$Cumul.eigen.values, lty = 2, cex = 1.5)
		}
		
		if(info == "area"){
			
			plot(x@PSR, 
					xlab = "cumulative eigenvalues(%)", 
					ylab = "r squared", ylim = c(0,1), 
					xlim = c(0,1), bty = "L", expected = segments(0, 0, 1, 1), pch = 16, ...)
		}
	}
}