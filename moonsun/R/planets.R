`planets` <-
function (jday=jd(),show.sun=TRUE,show.moon=TRUE) 
{
	res = rbind(
		mercury(jday),
		venus(jday),
		mars(jday),
		jupiter(jday),
		saturn(jday),
		uranus(jday),
		neptune(jday),
		pluto(jday));
	if (show.moon) res=rbind(moon(jday),res);
	if (show.sun) res=rbind(sun(jday),res);
	return(res);
}

