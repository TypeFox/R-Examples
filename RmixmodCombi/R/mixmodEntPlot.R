mixmodEntPlot <-
function(z, combiM, abc = c("standard", "normalized"), reg = c(2), ICL = 0, ...)
{
	oldpar <- par(no.readonly = TRUE)
	on.exit(par(oldpar))

	ent <- numeric()
	Kmax <- ncol(z)
	z0 <- z
	for (K in Kmax:1) 
		{
			z0 <- t(combiM[[K]] %*% t(z0))
			ent[K] <- -sum(xlog(z0))
		}

	if (any(abc == "normalized"))
		{
			mergedn <- numeric()
			z0 <- z
			for (K in (Kmax-1):1)
				{
					z0 <- t(combiM[[K+1]] %*% t(z0))
					mergedn[K] = sum(sapply(mixmodMap_M2V(mixmodMap(z0)),function(x) any(which(as.logical(combiM[[K]][rowSums(combiM[[K]])==2,]))==x)))				}
		}
		
	if (Kmax == 2) reg = NULL
		
	if (any(abc == "standard"))
		{
			par(oma=c(0,0,4,0), mar = c(5,4,0,2)+0.1)
			plot(1:Kmax, ent, xlab = "Number of clusters", ylab = "Entropy", xaxp = c(1, Kmax, Kmax-1), xaxt = "n", ...)
			if (any(reg == 2)) 
				{	
					pcwsreg <- pcws2_reg(1:Kmax,ent)
					lines(1:pcwsreg$c, pcwsreg$a1*(1:pcwsreg$c) + pcwsreg$b1, lty = 2, col = "red")
					lines(pcwsreg$c:Kmax, pcwsreg$a2*(pcwsreg$c:Kmax) + pcwsreg$b2, lty = 2, col = "red")
				}
			if (any(reg == 3)) 
				{	
					pcwsreg <- pcws3_reg(1:Kmax,ent)
					lines(1:pcwsreg$c1, pcwsreg$a1*(1:pcwsreg$c1) + pcwsreg$b1, lty = 2, col = "blue")
					lines(pcwsreg$c1:pcwsreg$c2, pcwsreg$a2*(pcwsreg$c1:pcwsreg$c2) + pcwsreg$b2, lty = 2, col = "blue")
					lines(pcwsreg$c2:Kmax, pcwsreg$a3*(pcwsreg$c2:Kmax) + pcwsreg$b3, lty = 2, col = "blue")
				}
			if (ICL > 1) axis(1, at = 1:(ICL-1), labels = as.character(1:(ICL - 1)), col.ticks = "black")
			if (ICL > 0) axis(1, at = ICL, labels = "ICL", col.ticks = "red")
			if (ICL < Kmax - 1 & ICL > 0) axis(1, at = (ICL + 1):(Kmax - 1), labels = as.character((ICL + 1) : (Kmax - 1)), col.ticks = "black")
			if (ICL == 0) axis(1, at = 1:Kmax, labels = as.character(1:Kmax)) 
			
			title("Entropy plot", outer=TRUE)

			par(ask=TRUE)
		}
	if (any(abc == "normalized"))
		{
			par(oma=c(0,0,4,0), mar = c(5,4,0,2)+0.1)
			plot(cumsum(c(0,mergedn)), ent, xlab = "Cumul. count of merged obs.", ylab = "Entropy", ...)
			X <- cumsum(c(0,mergedn))
			if (any(reg == 2)) 
				{	
					pcwsreg <- pcws2_reg(X,ent)
					lines(X[1:pcwsreg$c], pcwsreg$a1*(X[1:pcwsreg$c]) + pcwsreg$b1, lty = 2, col = "red")
					lines(X[pcwsreg$c:Kmax], pcwsreg$a2*(X[pcwsreg$c:Kmax]) + pcwsreg$b2, lty = 2, col = "red")
				}
			if (any(reg == 3)) 
				{	
					pcwsreg <- pcws3_reg(X,ent)
					lines(X[1:pcwsreg$c1], pcwsreg$a1*(X[1:pcwsreg$c1]) + pcwsreg$b1, lty = 2, col = "blue")
					lines(X[pcwsreg$c1:pcwsreg$c2], pcwsreg$a2*(X[pcwsreg$c1:pcwsreg$c2]) + pcwsreg$b2, lty = 2, col = "blue")
					lines(X[pcwsreg$c2:Kmax], pcwsreg$a3*(X[pcwsreg$c2:Kmax]) + pcwsreg$b3, lty = 2, col = "blue")
					par
				}
			axis(1, at = X[ICL], labels = "ICL", col.ticks = "red", padj = 1.5)
			title("Normalized entropy plot", outer=TRUE)
		}

}
