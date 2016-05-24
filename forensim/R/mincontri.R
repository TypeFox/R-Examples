
#Hinda Haned,  december 10th 2008



mincontri<-
function(mix,loc=NULL)
{ 
if(is.null(loc))
{
	tmp <- mix@mix.all
	tmp2 <- max(sapply(tmp,length))/2
	return(ceiling(tmp2))
}

else{
	tmp <- mix@mix.all[loc]
	tmp2 <- max(sapply(tmp,length))/2
	return(ceiling(tmp2))
}

}


