`as.eqc` <-
function (x,time=lst(),phi=getOption("latitude")) 
{
	time = as.vector(time);
	r = 180/pi;

	if (inherits(x,"hoc")) {
			a = x$az/r; alpha = x$alt/r; phi = phi/r;
			
			delta = asin(sin(alpha)*sin(phi)+cos(alpha)*cos(phi)*cos(a));
			h = acos((sin(alpha)-sin(phi)*sin(delta))/(cos(phi)*cos(delta)));

			delta = delta*r;
			h = h*r;
			h[sin(a)>0] = 360 - h[sin(a)>0];
			h[h == 360] = 0;

			h = time - h/15;
			h[h<0] = h[h<0] + 24;
	
			res = data.frame(ra=as.vector(h),d=as.vector(delta));
			rownames(res)=rownames(x);
			class(res$ra)="time";
			class(res$d)="dms";
			class(res)=c("eqc","apos","data.frame");
			return(res);
		}
	else if (inherits(x,"ecc")) {
	 		r = 180/pi;
			l = x$lat/r; b=x$long/r; 

        		a = atan2(sin(l)*0.917464059943984-tan(b)*0.3978186756690804,cos(l));
       		d = asin(sin(b)*0.917464059943984+cos(b)*0.3978186756690804*sin(l));

        		a = a*r/15;
			a = a %% 24;

        		d = d*r;

        		res=data.frame(ra=as.vector(a),d=as.vector(d));
        		rownames(res)=rownames(x);
			class(res$ra)="time";
			class(res$d)="dms";
        		class(res)=c("eqc","apos","data.frame");
        		return(res);
		}
	else {
			res = as.data.frame(x);
			names(res) = c("ra","d");
			class(res)=c("eqc","apos","data.frame");
			class(res$ra)="time";
			class(res$d)="dms";
			return(res);
	}
}

