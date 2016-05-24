print.lars <-
function(x, ...)
{
	cat("\nCall:\n")
	dput(x$call)
	cat("R-squared:", format(round(rev(x$R2)[1], 3)), "\n")
	actions <- x$actions
	jactions <- unlist(actions)
	jsteps <- rep(seq(along = actions), sapply(actions, length))
	actmat <- rbind(jsteps, jactions)
	vn <- names(jactions)
	if(is.null(vn))
		vn <- rep("", length(jactions))
	dimnames(actmat) <- list(c("Step", "Var"), vn)
	cat(paste("Sequence of", x$type, "moves:\n"))
	print(actmat[2:1,  ])
	invisible(x)
}

