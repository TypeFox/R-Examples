unif.test.statistic <- function(x, type, transform=NULL) {


if (!is.null(transform)) {
	if (transform=="spacings") {
		ordered.sample <- sort(x)
		ordered.sample <- c(0, ordered.sample, 1)
		spacings <- diff(ordered.sample)
		ordered.spacings <- sort(spacings)
		n <- length(x)
		S0 <- (n+1)*ordered.spacings[1]
		S <- (n+1-(1:n))*diff(ordered.spacings)
		V <- cumsum(c(S0, S))
		x <- V
	}
}

if (type=="greenwood") {
	
	ordered.sample <- sort(x)
	ordered.sample <- c(0, ordered.sample, 1)
	spacings <- diff(ordered.sample)
	stat <- sum(spacings^2)
	
} else if (type=="spacings.max") {
	
	ordered.sample <- sort(x)
	ordered.sample <- c(0, ordered.sample, 1)
	spacings <- diff(ordered.sample)
	n <- length(x)
	stat <- max(abs(spacings-1/n))
}

else if (type=="qm") {    # Quesenberry-Miller

	ordered.sample <- sort(x)
	n.sample <- length(x)
	ordered.sample <- c(0, ordered.sample, 1)
	spacings <- diff(ordered.sample)
	stat <- sum(spacings^2)+sum(spacings[1:n.sample]*spacings[2:(n.sample+1)])

} else if (type=="ks") {

	ordered.sample <- sort(x)
	n.sample <- length(x)
	Dplus <- max((1:n.sample)/n.sample - ordered.sample)
	Dminus <- max(ordered.sample - (0:(n.sample-1))/n.sample)
	stat <- max(Dplus, Dminus)

} else if (type=="V") {

	ordered.sample <- sort(x)
	n.sample <- length(x)
	Dplus <- max((1:n.sample)/n.sample - ordered.sample)
	Dminus <- max(ordered.sample - (0:(n.sample-1))/n.sample)
	stat <- Dplus + Dminus

} else if (type=="cvm") {

	ordered.sample <- sort(x)
	n.sample <- length(x)
	stat <- sum((ordered.sample - (seq(from = 1, to = (2 * n.sample - 1), by = 2)/(2 * n.sample)))^2) + 1/(12 * n.sample)

} else stop("The goodness-of-fit test must be one of: Greenwood, Quesenberry-Miller, Kolmogorov-Smirnov, V=D+ +  D-, or Cramer-Von Mises")


return(stat)

}