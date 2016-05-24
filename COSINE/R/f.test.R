f.test <-
function(data,type)
{
	lg <- length(type)
	n <- tapply(type,type,length)
	g <- length(n)
	x <- tapply(data,type,mean)
	y <- tapply(data,type,var)
	xx <- mean(data)
	a <- sum(n*(x-xx)^2)/(g-1)
	b <- sum((n-1)*y)/(lg-g)
	s <- a/b
	return(s)
}

