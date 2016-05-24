summary.tosANYN <-
function (object, quiet = FALSE, ...) 
{
    hwtosop <- object
    ntests <- length(hwtosop$allpvals)
    if (quiet == FALSE) {
        cat("There are ", ntests, " hypothesis tests altogether\n")
        cat("There were ", hwtosop$nreject, " reject(s)\n")
	cat("P-val adjustment method was: ", object$mc.method[1], "\n")
    }

    rejix <- which(object$allpvals < object$alpha)

    vmat <- cbind(object$allbigscale[rejix]-1,
		object$alllitscale[rejix],
		object$alllv[rejix],
		object$allindex[rejix])

    v <- as.data.frame(t(vmat))
    class(v) <- "list" 


    if (quiet==FALSE)	{
	if (hwtosop$nreject != 0)	{
	  cat("Listing rejects...\n")
	  for(i in 1:hwtosop$nreject)	{
		cat("P: ", v[[i]][1], " HWTlev: ", v[[i]][2],
			" Max Poss Ix: ", v[[i]][3],
			" Indices: ", v[[i]][c(-1,-2,-3)], "\n") 
		}
	
	  }
	}

    vret <- list(rejlist = v, nreject = hwtosop$nreject, mctype = object$p.adjust.method[1])

    return(invisible(vret))
}
