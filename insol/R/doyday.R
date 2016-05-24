doyday <-
function(year,doy){
	if (nargs() < 1 ) {cat("USAGE: doyday(year,doy) \nUSAGE: doyday \n(year.dy)"); return()}
	fyy=floor(year)
	if (nargs() == 1 ) {
	nd=ifelse((fyy%%4==0 & fyy%%100!=0) | fyy%%400==0,366,365)
	dd=(year-fyy)*nd
	} else {
	yy=year
	dd = doy }
	fdd=floor(dd)
	hh=(dd-fdd)*24
	fhh=floor(hh)
	mm=(hh-fhh)*60
	fmm=floor(mm)
	ss=(mm-fmm)*60
	fss=floor(ss)
	return(strptime(paste(fyy,fdd,fhh,fmm,fss),format="%Y %j %H %M %S"))
}

