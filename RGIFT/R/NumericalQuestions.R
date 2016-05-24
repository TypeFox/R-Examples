GIFTNQ<-function(qtxt, ans, err=0)
{

	if(!(length(ans)==1 | length(ans==2)))
		stop("Wrong number of answers.")
	cat(qtxt, "{", sep="")

	if(length(ans)==1)
		cat("#", ans, ":", err, "}.\n\n", sep="")

	if(length(ans)==2)
		cat("#", ans, ":", err[1], "..", err[2], "}.\n\n", sep="")

}

