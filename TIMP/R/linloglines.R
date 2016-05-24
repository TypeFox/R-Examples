"linloglines" <-
function (x, mu, alpha) 
{
    if(alpha != 0) {
     for (i in 1:length(x)) {
        if ((x[i] - mu)/alpha >= 1) 
            x[i] <- alpha + (alpha * log10((x[i] - mu)/alpha))
        else if ((x[i] - mu)/alpha < -1) 
            x[i] <- -alpha - (alpha * log10((mu - x[i])/alpha))
        else x[i] <- x[i] - mu
     }
    }
    else    
	     for (i in 1:length(x))
		  x[i] <- x[i] - mu
    x
}

