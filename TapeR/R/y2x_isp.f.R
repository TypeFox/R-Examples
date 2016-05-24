y2x_isp.f <-
function(x, x.grd, y.grd, ...){
#   ------------------------------------------------------------------------------------------------

		 y = spline(xout = c(x), x = x.grd, y = y.grd, method = "natural", ties = mean)$y
 	#	 y = predict(interpSpline( x = x.grd, y = y.grd),x=x)$y

		 return(y*y)
	}
