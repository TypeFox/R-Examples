`sun` <-
function (jday=jd()) 
{
	d = as.vector(jday) - 2447891.5;
	n = 360/365.242191*d
	n = n %% 360;
	m = n - 3.365119;
	m = m %% 360;
	ec = 360/pi*0.016713*sin(m/180*pi);
	lambda = n + ec + 279.403303;
	lambda = lambda %% 360;

	v = m + ec

	r = 0.9997207/(1+0.016713*cos(v/180*pi))
	th = 0.533128/r

	nam = c();
	for (i in 1:length(jday))
		{ nam = c(nam,paste(format.jd(jday[i]),"Sun",sep="-"));
		}

	pos = ecc(lambda,0,nam);
	pos = as.eqc(pos);
	pos$phase=NA
	pos$angle=NA
	pos$dist=round(r,2);
	class(pos$dist)="numeric"
	pos$size=round(th,1);
	class(pos$size)="numeric"
	pos$mag=NA;
	
	return(pos);
}

