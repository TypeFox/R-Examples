#Converts a matrix of points at once
angles <- function(vec)
{
	d <- ncol(vec)
	cart <- matrix(0,nrow=nrow(vec),ncol=ncol(vec))
	cosine <- cos(vec[,-1])
	sine   <- sin(vec[,-1])
	if(d == 2)
	{
		cart[,1] <- vec[,1] * cosine
		cart[,2] <- vec[,1] * sine
	}
	else if(d == 3)
	{
		cart[,1] <- vec[,1] * cosine[,1]
		cart[,2] <- vec[,1] * sine[,1] * cosine[,2]
		cart[,3] <- vec[,1] * sine[,1] * sine[,2]
	}
	else
	{
		cart[,1] <- vec[,1] * cosine[,1]
		cart[,2] <- vec[,1] * sine[,1] * cosine[,2]
		for(i in 3:(d-1))	#For d coordinates
		{
			cart[,i] <- vec[,1] * apply(sine[,1:(i-1)],FUN="prod",MARGIN=1) * cosine[,i]
		}
		cart[,d] <- vec[,1] * apply(sine[,1:(d-1)],FUN="prod",MARGIN=1)
	}
	return(cart)
}
