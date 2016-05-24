GIFTMW<-function(qtxt1, qtxt2, anstxt, rightans)
{
	cat(qtxt1, "{\n", sep="")

	for(i in 1:length(anstxt))
	{

		if(i==rightans)
			cat("=", anstxt[i], "\n", sep="")
		else
			cat("~", anstxt[i], "\n", sep="")

	}
	cat("}\n",qtxt2, "\n\n")
}
