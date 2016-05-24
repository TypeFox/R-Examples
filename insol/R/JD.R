JD <-
function(x, inverse=FALSE) {
# 	options(digits=12)
	if (inverse){
		return(as.POSIXct((x-2440587.5)*86400,origin=ISOdate(1970,01,01,0,0,0),format="%Y-%m-%d %H:%M:%S" ))
	}else{
		return(as.numeric(x)/86400 + 2440587.5)
	}
}

