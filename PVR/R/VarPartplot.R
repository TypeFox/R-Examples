VarPartplot <- function(x, ...){
	
	A <- round(x@VarPart[[1]], 5)
	B <- round(x@VarPart[[2]], 5)
	C <- round(x@VarPart[[3]], 5)
	D <- round(x@VarPart[[4]], 5)
	
	plot(c(0.5,0.5) ~ c(0,1), xlim = c(-0.01, 1.01), ylim = c(-0.01, 1.01), 
			xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "", 
			type = "l", lwd = 2)
	segments(0, 0.5, 0, 0.7, lty = 2); 
	segments(A, 0.5, A, 0.7, lty = 2); 
	segments((A + B), 0.5, (A + B), 0.7, lty = 2); 
	segments((A + B + C), 0.5, (A + B + C), 0.7, lty = 2)
	coordsAB <- data.frame(c(0, (A+B), (A+B), 0), c(0.625,  0.625, 0.5, 0.5))
	polygon(coordsAB, col = "white")
	text((A/2), 0.6625, "a")
	coordsBC <- data.frame(c(A, (A+B+C), (A+B+C), A), c(0.375,  0.375, 0.5, 0.5))
	text((A+(B/2)), 0.6625, "b")
	polygon(coordsBC, col = "white")
	text((A+B+(C/2)), 0.6625, "c")
	text((A+B+C+(D/2)), 0.6625, "d")
	lines(c(0.5,0.5) ~ c(0,1), lwd = 2)
	
	cat("[a] - ", "Variation explained by environmental variables", "\n", 
			"[b] - ", "Shared variation between environmental variables and phylogeny(PVR)", "\n",
			"[c] - ", "Variation explained by phylogeny(PVR)", "\n",
			"[d]     - ", "Unexplained variation", "\n", sep = "")
}

