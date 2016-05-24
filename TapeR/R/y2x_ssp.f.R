y2x_ssp.f <-
function(x,x.grd,y.grd, ...){ssp = smooth.spline(x.grd,y.grd);return(predict(ssp,x)$y^2)}
