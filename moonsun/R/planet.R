`planet` <-
function (jday=jd(),name="",inner=FALSE,tp,ep,oo,e,a,i,om,th,mag) 
{
	d = as.vector(jday) - 2447891.5;

	np = 360/365.242191*d/tp;
	np = np %% 360;
	mp = np + ep - oo;
	l = np + 360/pi*e*sin(mp/180*pi)+ep;
	l = l %% 360;
	vp = l - oo;
	r = a*(1-e*e)/(1+e*cos(vp/180*pi));

	ne = 360/365.242191*d/1.00004;
	ne = ne %% 360;
	me = ne + 99.403308 - 102.768413;
	L = ne + 360/pi*0.016713*sin(me/180*pi)+99.403308;
	L = L%% 360;
	ve = L - 102.768413;
	R = 0.9997206756/(1+0.016713*cos(ve/180*pi));

	psi = asin(sin((l-om)/180*pi)*sin(i/180*pi))*180/pi;
	y = sin((l-om)/180*pi)*cos(i/180*pi);
	x = cos((l-om)/180*pi);
	
	lp = atan2(y,x)*180/pi+om;
	rp = r*cos(psi/180*pi);

	if (inner) {
		A = atan2(rp*sin((L-lp)/180*pi),R-rp*cos((L-lp)/180*pi))*180/pi;
		lambda = 180 + L + A;
		}
	else lambda = atan2(R*sin((lp-L)/180*pi),rp-R*cos((lp-L)/180*pi))*180/pi+lp;
 

	lambda = lambda %% 360;
	beta = atan2(rp*tan(psi/180*pi)*sin((lambda-lp)/180*pi),R*sin((lp-L)/180*pi))*180/pi;
	beta[beta< -90] = beta[beta< -90] + 180;
	beta[beta>90] = beta[beta>90] - 180;
	nam = c();
        for (i in 1:length(jday))
                { nam = c(nam,paste(format.jd(jday[i]),name,sep="-"));
                }
	
	pos = ecc(lambda,beta,nam);
	pos = as.eqc(pos);
	pos$phase = 100*round(0.5*(1+cos((lambda-l)/180*pi)),3);

	spos = sun(jday);

	pos$angle = round(180/pi*atan2(cos(spos$d/180*pi)*sin(15*(spos$ra-pos$ra)/180*pi),
		cos(pos$d/180*pi)*sin(spos$d/180*pi)-sin(pos$d/180*pi)*cos(spos$d/180*pi)*cos(15*(spos$ra-pos$ra)/180*pi)),1)	

	pos$dist = round(sqrt(R*R+r*r-2*r*R*cos((l-L)/180*pi)),2)

	pos$size = round(th/pos$dist,1);
	pos$mag = round(5*log10(r*pos$dist/sqrt(pos$phase/100))+mag,1);

	class(pos$angle)="numeric"
	class(pos$dist)="numeric"
	class(pos$mag)="numeric"

	return(pos);
}

