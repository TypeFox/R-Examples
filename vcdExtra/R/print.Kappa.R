# Print method for Kappa: Add a column showing z values 
## DONE: now set digits
## DONE: now include CI  
print.Kappa <-
		function (x, digits=max(getOption("digits") - 3, 3), CI=FALSE, level=0.95, ...) 
{
	tab <- rbind(x$Unweighted, x$Weighted)
	z <- tab[,1] / tab[,2]
	tab <- cbind(tab, z)
	if (CI) {
		q <-  qnorm((1 + level)/2)
		lower <- tab[,1] - q * tab[,2]
		upper <- tab[,1] + q * tab[,2]
		tab <- cbind(tab, lower, upper)
	}
	
	rownames(tab) <- names(x)[1:2]
	print(tab, digits=digits, ...)
	invisible(x)
}
