library('mapmisc')
Spackages = c('rgdal', 'rgeos', 'geosphere', 'maptools')

if(all(unlist(mapply(requireNamespace, package=Spackages, MoreArgs=list(quietly=TRUE))))){
	library('rgdal')

	data('wrld_simpl', package='maptools')
	
	country='Japan'
	Dcountry  = grep(country, wrld_simpl$NAME)
	x=wrld_simpl[Dcountry,]

	myCrsO = moll(x, angle=25)
	
	xTcrop = wrapPoly(x=wrld_simpl, crs=myCrsO)
	DcountryT  = grep(country, xTcrop$NAME)
	
	map.new(xTcrop)
	plot(attributes(myCrsO)$ellipse, add=TRUE, col='lightBlue', border='blue')
	plot(xTcrop,add=TRUE)
	plot(xTcrop[DcountryT,], col='red', add=TRUE)
	
	gridlinesWrap(myCrsO, lty=2, col='orange')
	
	
	
	country='Madagascar'
	Dcountry  = grep(country, wrld_simpl$NAME)
	x=wrld_simpl[Dcountry,]
	
	myCrsMoll = moll(x, angle=100)
	
	
	plot(wrld_simpl)
	plot(attributes(myCrsMoll)$crop, col='red', border='yellow', add=TRUE)
	
	xTcrop = wrapPoly(x=wrld_simpl, crs=myCrsMoll)
	DcountryT  = grep(country, xTcrop$NAME)
	
	map.new(xTcrop)
	plot(attributes(myCrsO)$ellipse, add=TRUE, col='lightBlue', border='blue')
	plot(xTcrop,add=TRUE)
	plot(xTcrop[DcountryT,], col='green', add=TRUE)
	
	gridlinesWrap(crs=xTcrop, lty=2, col='red')

}