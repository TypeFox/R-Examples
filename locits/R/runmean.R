runmean <-
function(x, binwidth){

lx <- length(x)

if (binwidth != round(binwidth))
	stop("binwidth should be an integer")
else if (binwidth < 1)
	stop("binwidth needs to be positive")

if (binwidth==1)
	return(x)
else if (binwidth==2)
	xnew <- c(x[1], x)
else if (binwidth==3)
	xnew <- c(x[1], x, x[length(x)])

else if (binwidth %% 2 == 1)	{
	bw2 <- (binwidth-1)/2
	xnew <- c(x[(bw2-1):1], x, x[lx:(lx-bw2)])
	}
else 	{
	bw2 <- binwidth/2
	xnew <- c(x[(bw2-1):1], x, x[lx:(lx-bw2+1)])
	}

running.mean(xnew, binwidth=binwidth)
}
