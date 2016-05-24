`dprime.ABX` <-function(Hits, FA, zdiff, Pc.unb, 
				method = "diff") {
	Call <- match.call() 
	 if (pmatch("Hits", names(Call), 0) > 0) {
        if (pmatch("FA", names(Call), 0) > 0) {
                        zdiff <- qnorm(Hits) - qnorm(FA)
                        Pc.unb <- pnorm(zdiff/2)
                } else {
                	 zdiff <- qnorm(Hits) - qnorm(1-Hits)
                     Pc.unb <- pnorm(zdiff/2)
        } } else {
        	if (pmatch("zdiff", names(Call), 0) > 0)
                   { Pc.unb <- pnorm(zdiff/2)       } 
                        }
	if (Pc.unb < 0.5) stop("Only valid for Pc.unb > 0.5")
	root2 <- sqrt(2)
	if (method == "diff") {
		root6 <- sqrt(6)
		est.dp <- function(dp) {
			Pc.unb - pnorm(dp/root2)*pnorm(dp/root6) - 
				 pnorm(-dp/root2)*pnorm(-dp/root6)
		}
     	dp.res <- uniroot(est.dp, interval=c(0,10))
     	dprime <- dp.res$root             
	 } else 
	 {
	 if (method == "IO") {
	 	est.dp <- function(dp) {
		Pc.unb - pnorm(dp/root2)*pnorm(dp/2) - 
			 pnorm(-dp/root2)*pnorm(-dp/2)
		}
     dp.res <- uniroot(est.dp, interval=c(0,10))
     dprime <- dp.res$root             
	} else
	{stop("method not defined; must be either diff or IO") } }
dprime
}
