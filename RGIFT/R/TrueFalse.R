GIFTTF<-function(qtxt, ans)
{
	if(!is.logical(ans))
		stop("Answer is not logical.")

        cat("\n",qtxt, "{", as.character(ans), "}\n", sep="")
}

