
GIFTMC<-function(qtxt, anstxt, rightans=1, wright=NULL, wwrong=NULL)
{

	if(is.null(wright))
		wright<-"100"

	if(is.null(wwrong))
		wwrong<-"0"


	#Split weight of right answer among all right answers
	if(length(rightans)>1)
	{
		if(length(rightans)!=length(wright))
			stop("Number of right answers and weights differ.")
	}

        cat("\n",qtxt, "\n", sep="")

        cat("{\n", sep="")

        for(i in 1:length(anstxt))
        {
		if(i %in% rightans)
		{
			idx<-which(i==rightans)
			if(wright[idx]=="100")#if(idx==1)
	                	cat("=", anstxt[i],"\n", sep="")
			else
		                cat("~%", wright[idx],"%", anstxt[i],"\n", sep="")
		}
		else
		{
	                cat("~%", wwrong,"%", anstxt[i],"\n", sep="")
		}
        }

        cat("}\n", sep="")

}
