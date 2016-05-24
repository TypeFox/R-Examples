GIFTSA<-function(qtxt, anstxt, wright="100")
{


        #Split weight of right answer among all right answers
        if(length(anstxt)>1)
	{
		if(length(anstxt)!=length(wright))
		warning("Number of right answers and weights differ.")
		wright<-rep(wright, length.out=length(anstxt))
	}
	else
	{
		wright<-rep(wright, length(anstxt))
	}



	cat(qtxt," {\n", sep="")

	for(i in 1:length(anstxt))
	{
#		if(wright[i]=="100")
#			cat("=", anstxt[i], "\n", sep="")
#		else
			cat("=%", wright[i], "%", anstxt[i], "\n", sep="")
	}

	cat("}\n\n")


}
