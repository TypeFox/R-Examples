`dprime.SD` <-
function (H, FA, zdiff, Pcmax, method = "diff") {
	if (method == "diff") {
		if (isTRUE(all.equal(H, FA))) return(0)
		root2 <- sqrt(2)
		k <- root2 * qnorm(FA/2)
		est.dp <- function(dp)
			{ H - pnorm((k + dp)/root2) - pnorm((k - dp)/root2) }
		dp.res <- uniroot(est.dp, interval = c(0, 10))
		dprime <- dp.res$root
	} else
	if (method == "IO")	{
		Call <- match.call() 
    if (pmatch("H", names(Call), 0) > 0) {
    	if (pmatch("FA", names(Call), 0) > 0) {
			zdiff <- qnorm(H) - qnorm(FA)
			Pcmax <- pnorm(zdiff/2)
		} else {
			zdiff <- qnorm(H) - qnorm(1-H)
			Pcmax <- pnorm(zdiff/2)
	} } else {
			if (pmatch("zdiff", names(Call), 0) > 0)
				{ Pcmax <- pnorm(zdiff/2)	} 
			} 
	dprime <- sign(Pcmax - 0.5) * if ( Pcmax < 0.5 ) 2 * qnorm(0.5 * (1 + sqrt(2 * (1 - Pcmax) - 1))) else 2 * qnorm(0.5 * (1 + sqrt(2 * Pcmax - 1)))
#2 * qnorm(0.5 * (1 + sqrt(2*Pcmax - 1)))
	} else
	{stop("method must be one of diff or IO") }
	if(dprime < 0) warning("FA > H giving d' < 0!")
	dprime
}
