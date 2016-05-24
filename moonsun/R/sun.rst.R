`sun.rst` <-
function (jday=jd(),  phi = getOption("latitude"))
{
	
	rise = 0 
	for (i in 1:5) {
		rise[!is.finite(rise)] = 0
            rise = as.gmt(rst(sun(jday+rise/24),phi)$rise)
		}

	transit = 0 
	for (i in 1:5) {
		transit[!is.finite(transit)] = 0
		transit = as.gmt(rst(sun(jday+transit/24),phi)$transit)
		}

	set = 0 
	for (i in 1:5)  {
		set[!is.finite(set)] = 0
		set = as.gmt(rst(sun(jday+set/24),phi)$set)
		}
			
      res = data.frame(rise=as.vector(as.lst(rise)),
				transit=as.vector(as.lst(transit)),
			      set=as.vector(as.lst(set)))

      rownames(res) = format.jd(jday);

      class(res) = c("rst","data.frame");
      class(res$rise) = c("lst","time");
      class(res$set) = c("lst","time");
      class(res$transit) = c("lst","time");

      return(res);

	
}

