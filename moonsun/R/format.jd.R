`format.jd` <-
function (x,  ...)
{
	jd = x
	res = c()

for (ii in 1:length(jd)) 
 { 

	jj = jd[ii] + 0.5;
	i = floor(jj);
	f = jj - i;
	
	if (i>2299160) {
		a = floor((i -1867216.25)/36524.25);
		b = i + 1 + a - floor(a/4);
			}
	else b = i;

	c = b + 1524;
	d = floor((c-122.1)/365.25);
	e = floor(365.25*d);
	g = floor((c-e)/30.6001);
	day = c - e + f - floor(30.6001*g);

	if (g<13.5) month = g - 1
	else month = g - 13;

	if (month>2.5) year = d - 4716
	else year = d - 4715;

	jj = jj - 0.5;
	res=c(res,paste(year,"-",formatC(month,width=2,flag="0"),"-",formatC(floor(day),width=2,flag="0"),sep=""));
 }
	return(res);
}

