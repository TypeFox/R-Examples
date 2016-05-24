"dtt" <-
function (x, type=c("dct","dst","dht"), variant=2, inverted=FALSE) 
{
    if (is.matrix(x)) 
        t(apply(x, 1, dtt, type=type, variant=variant, inverted=inverted))
    else if (is.data.frame(x)) 
        t(apply(x, 1, dtt, type=type, variant=variant, inverted=inverted))
    else

{

	n <- length(x);
	res <- x;
	type <- match.arg(type);
	transform <- switch(type,dct=1,dst=2,dht=3);

if (transform == 1) {

	if ((variant == 1)) {
	for (k in 0:(n-1)) {
		res[k+1] <- 0.5*(x[1]+((-1)^k)*x[n]) + sum(x[2:(n-1)]*cos(pi/(n-1)*(1:(n-2))*k));
	}
	if (inverted) { res=res*(2/(n-1)); }
	}

	if ((variant == 2 && !inverted) || (variant == 3 && inverted)) {
	for (k in 0:(n-1)) {
		res[k+1] <- sum(x*cos(pi/n*((0:(n-1))+0.5)*k));
	}
	if (inverted) { res=res*(2/n); }
	}

	if ((variant == 3 && !inverted) || (variant == 2 && inverted)) {
	for (k in 0:(n-1)) {
		res[k+1] <- 0.5*x[1] + sum(x[2:n]*cos(pi/n*(1:(n-1))*(k+0.5)));
	}
	if (inverted) { res=res*(2/n); }
	}

	if ((variant == 4)) {
	for (k in 0:(n-1)) {
		res[k+1] <- sum(x*cos(pi/n*((0:(n-1))+0.5)*(k+0.5)));
	}
	if (inverted) { res=res*(2/n); }
	}
}
else if (transform == 2) {

	if (variant == 1) {
	for (k in 0:(n-1)) {
		res[k+1] <- sum(x*sin(pi/(n+1)*((0:(n-1))+1)*(k+1)));
	}
	if (inverted) { res=res*(2/(n-1)); }
	}

	if ((variant == 2 && !inverted) || (variant == 3 && inverted)) {
	for (k in 0:(n-1)) {
		res[k+1] <- sum(x*sin(pi/n*((0:(n-1))+0.5)*(k+1)));
	}
	if (inverted) { res=res*(2/n); }
	}

	if ((variant == 3 && !inverted) || (variant == 2 && inverted)) {
	for (k in 0:(n-1)) {
		res[k+1] <- (-1)^k/2*x[n] + sum(x[1:(n-1)]*sin(pi/n*((0:(n-2))+1)*(k+0.5)));
	}
	if (inverted) { res=res*(2/n); }
	}

	if (variant == 4) {
	for (k in 0:(n-1)) {
		res[k+1] <- sum(x*sin(pi/n*((0:(n-1))+0.5)*(k+0.5)));
	}
	if (inverted) { res=res*(2/n); }
	}

}

else if (transform == 3) {
        for (k in 0:(n-1)) {
                res[k+1] <- sum(x*(cos(2*pi/n*(0:(n-1))*k)+sin(2*pi/n*(0:(n-1))*k)));
        }
        if (inverted) { res=res*(1/n); }
}

	return(res);
}

}

