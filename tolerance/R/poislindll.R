poislind.ll <- function(x, theta = NULL, ...){
    if (class(x) != "table") {
        x <- table(x)
    }
	N <- as.numeric(x)
	x <- as.numeric(names(x))
	x.bar <- sum(N*x)/sum(N)
	if(is.null(theta)) theta <- (-(x.bar-1)+sqrt((x.bar-1)^2+8*x.bar))/(2*x.bar)
	ll.f <- function(theta) -2*sum(N)*log(theta)-sum(N*(log(x+theta+2)-log(theta+1)*(x+3)))
	fit <- try(suppressWarnings(stats4::mle(ll.f,start=list(theta=theta),lower=0,...)),silent=TRUE)
	if(class(fit)=="try-error"){
		fit <- try(suppressWarnings(stats4::mle(ll.f,start=list(theta=.8*theta),lower=0,...)),silent=TRUE)
		if(class(fit)=="try-error") stop(paste("Difficulty optimizing the MLE -- must try a different starting value for theta.","\n"))
	}
	fit
}
