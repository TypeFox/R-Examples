GIFTM<-function(qtxt, group1, group2)
{
	if(length(group1)!=length(group2))
		stop("Different number of words to match.")


	cat(qtxt, "{\n", sep="")

	for(i in 1:length(group1))
	{
		cat("=",group1[i], " -> ", group2[i], sep="")
	}

	cat("}\n\n")

}
