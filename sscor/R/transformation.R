transformationcor <- function(s11,s12) {
	d <- 0.5 + sqrt((s11-0.5)^2 + s12^2)
	b <- d - s11
	a <- (2*d -1 )/(d*(1-d))
	r <- (a * s12 * b)/sqrt((s12^2 + b^2)^2 + (s12 * a * b)^2 ) 
	return(r)
}

trafofisher <- function(x) return(sign(x)*(1/sqrt(2)*asin((3*(1-sqrt(1-x^2))-2)/(sqrt(1-x^2)+1))+pi/2^(3/2)))

trafofisherinv <- function(y) return(sign(y)*2^(3/2)*sqrt(1-cos(sqrt(2)*y))/(3-cos(sqrt(2)*y)))
