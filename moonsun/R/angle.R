`angle` <-
function (x,y=NULL) 
{
	r = 180/pi;
	if (inherits(x,"eqc")) a = x[,1]*15/r else a = x[,1]/r;
	d = x[,2]/r;

	res = c();

	n = nrow(x);

	if (!is.null(y)) {
		if (inherits(y,"eqc")) a2 = y[,1]*15/r else a2 = y[,1]/r;
		d2 = y[,2]/r;

	for (i in 1:n) 
		res = c(res, acos(sin(d[i])*sin(d2[i])+cos(d[i])*cos(d2[i])*cos(a[i]-a2[i])))
		res = round(res*r,1)
	}
	else {
	
		for (i in 1:(n-1))
			for (j in (i+1):n)
				res = c(res,acos(sin(d[i])*sin(d[j])+cos(d[i])*cos(d[j])*cos(a[i]-a[j])));

		res = round(res*r,1);

		class(res)=c("angle","dist");
		attr(res, "Size") = n;
		attr(res, "Labels") = rownames(x);
		attr(res, "Diag") = FALSE;
		attr(res, "Upper") = FALSE;

	}	
	
		return(res);


}

