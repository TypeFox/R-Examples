`as.hoc` <-
function (x,time=lst(),phi=getOption("latitude")) 
{
	time = as.vector(time);

        if (is.null(phi)) 
                        {
                        phi = 0;
                        warning("Your latitude is not set in the environment, assuming it is equal to 0");
                        }

	if (inherits(x,"eqc")) {

		r = 180/pi;

		h = (time - x$ra)*15/r;
		phi=phi/r; 		
		d = x$d/r;

		alpha = asin(sin(d)*sin(phi)+cos(d)*cos(phi)*cos(h));
		a = acos((sin(d)-sin(phi)*sin(alpha))/(cos(phi)*cos(alpha)));

		alpha = alpha*r;
		a = a*r;

		a[sin(h)>=0] = 360 - a[sin(h)>=0];
		a[a == 360] = 0;

		res=data.frame(az=as.vector(a),alt=as.vector(alpha));
		rownames(res)=rownames(x);
		class(res$az)="dms";
		class(res$alt)="dms";
		class(res)=c("hoc","apos","data.frame");
		return(res);
	}
	else if (inherits(x,"ecc")) return(as.hoc(as.eqc(x)))
	else {
		res = as.data.frame(x);
		names(res) = c("az","alt");
		class(res$az)="dms";
		class(res$alt)="dms";
		class(res)=c("hoc","apos","data.frame");
		return(res);
		}

}

