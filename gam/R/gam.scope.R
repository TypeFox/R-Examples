"gam.scope" <-
function(frame, response = 1, smoother = "s", arg = NULL, form = TRUE)
{
	vnames <- names(frame)
	vnames <- vnames[ - response]
	step.list <- as.list(vnames)
	names(step.list) <- vnames
	for(vname in vnames) {
		junk <- c("1", vname)
		if(is.vector(frame[[vname]]))
			junk <- c(junk, paste(smoother, "(", vname, if(is.null(
				arg)) ")" else paste(",", arg, ")", sep = ""),
				sep = ""))
		if(form)
			junk <- eval(parse(text = paste("~", paste(junk, 
				collapse = "+"))))
		step.list[[vname]] <- junk
	}
	step.list
}
