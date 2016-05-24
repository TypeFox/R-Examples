`moon` <-
function (jday=jd()+gmt()/24) 
{
	d = as.vector(jday) - 2447891.5;
	n = 360/365.242191*d
	n = n %% 360;
	m = n - 3.365119;
	m = m %% 360;
	ec = 360/pi*0.016713*sin(m/180*pi);
	lambda = n + ec + 279.403303;
	lambda = lambda %% 360;

	l = 13.1763966*d+318.351648;
	l = l %% 360;
	mm = l-0.1114041*d-36.340410;
	mm = mm %% 360;
	n = 318.510107-0.0529539*d;
	n = n %% 360;
	ev = 1.2739*sin((2*(l-lambda)-mm)/180*pi);
	ae = 0.1858*sin(m/180*pi);
	a3 = 0.37**sin(m/180*pi);
	mpm = mm+ev-ae-a3;
	ec = 6.2886*sin(mpm/180*pi);
	a4 = 0.214*sin(2*mpm/180*pi);
	lp = l+ev+ec-ae+a4;
	v = 0.6583*sin(2*(lp-lambda)/180*pi);
	lpp = lp+v;
	np = n - 0.16*sin(m/180*pi);
	y = sin((lpp-np)/180*pi)*0.995970321;
	x = cos((lpp-np)/180*pi);
	lm = atan2(y,x)*180/pi+np;

	beta = asin(sin((lpp-np)/180*pi)*0.08968344185);

	nam = c();
	for (i in 1:length(jday))
		{ nam = c(nam,paste(format.jd(jday[i]),"Moon",sep="-"));
		}

	pos = ecc(lm,beta,nam);
	pos = as.eqc(pos);
	pos$phase = 100*round(0.5*(1-cos((lambda-lpp)/180*pi)),3);

      spos = sun(jday);

      pos$angle = round(180/pi*atan2(cos(spos$d/180*pi)*sin(15*(spos$ra-pos$ra)/180*pi),
                cos(pos$d/180*pi)*sin(spos$d/180*pi)-sin(pos$d/180*pi)*cos(spos$d/180*pi)*cos(15*(spos$ra-pos$ra)/180*pi)),1)    

	class(pos$angle)="numeric"

	pos$dist = round(0.996986/(1+0.0549*cos((mpm+ec)/180*pi)),2)
	class(pos$dist)="numeric"
	pos$size = round(60*0.5181/pos$dist,1)
	class(pos$size)="numeric"
	pos$mag = NA;

	return(pos);
}

