"wofz" <-
function(z){

ans <- .C("IPwofzvec",
		x=as.double(Re(z)),
		y=as.double(Im(z)),
		ansx=as.double(Re(z)),
		ansy=as.double(Im(z)),
		n = as.integer(length(z)),
		error = as.integer(0),
		PACKAGE="NORMT3")
if (ans$error != 0)
	stop(paste("Error code from TOMS 680 was ", ans$error))

return(complex(real=ans$ansx, imaginary=ans$ansy))

}
