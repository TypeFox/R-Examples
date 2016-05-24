`as.ecc` <-
function (x) 
{
	if (inherits(x,"eqc")) {
		r = 180/pi;
		a = x$ra*15/r; d=x$d/r; 

		l = atan2(sin(a)*0.917464059943984+tan(d)*0.3978186756690804,cos(a));
		b = asin(sin(d)*0.917464059943984-cos(d)*0.3978186756690804*sin(a));

		l = l*r;
		l = l %% 360;

		b = b*r;

		res=data.frame(lat=as.vector(l),long=as.vector(b));
		rownames(res)=rownames(x);
		class(res$lat)="dms";
		class(res$long)="dms";
		class(res)=c("ecc","apos","data.frame");
		return(res);
	} 
	else if (inherits(x,"hoc")) return(as.ecc(as.eqc(x)))
	else {
		res = as.data.frame(x);
		names(res)=c("lat","long");		
		class(res)=c("ecc","apos","data.frame");
		class(res$lat)="dms";
		class(res$long)="dms";
		return(res);
	}
}

