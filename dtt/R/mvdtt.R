"mvdtt" <-
function (x, type=c("dct","dst","dht"), variant=2, inverted=FALSE) 
{
    if (!is.matrix(x)) stop ("Parameter must be a matrix");

	s <- dtt(x, type=type, variant=variant, inverted=inverted)
	r <- dtt(t(s), type=type, variant=variant, inverted=inverted)
	res <- t(r);
	return(res);

}

