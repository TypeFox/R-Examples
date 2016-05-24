`jd` <-
function (year = NULL, month = NULL, day= NULL, epoch = Sys.time(), length = 1, by = 1 ) 
{
	if (any(is.null(year),is.null(month),is.null(day))) {
				daylt = as.POSIXlt(epoch);
				year = daylt$year+1900;
				month = daylt$mon+1;
				day = daylt$mday;
			}
			
	if ((month == 1) || (month == 2)) {
				year = year - 1;
				month = month + 12;
				}

	if ((year > 1582) ||
	     ((year == 1582) && (month > 10)) ||
	     ((year == 1582) && (month == 10) && (day >= 15))) {
			a = floor(year/100);
			b = 2 - a + floor(a/4);
	     }

	if (year<0) c = floor((365.25 * year) - 0.75)
		else 	c = floor (365.25 * year);

	d = floor(30.6001 * (month + 1));

	
	julian = seq(b+c+d+day+1720994.5,length = length,by=by);
	class(julian)="jd";
	
	return(julian);
}

