sphdist = function( long1, lat1, long2, lat2, degrees=FALSE) {
 	cf = 1.0
	if (degrees) cf = pi/180

	#--- Convert both points to rectangular coordinates. ---
        polrec = function(r,theta) r*c(cos(theta), sin(theta))

	tmp = polrec(1.0, lat1/cf);  rxy=tmp[1]; z1=tmp[2]
	tmp = polrec(rxy, long1/cf); x1 = tmp[1]; y1=tmp[2]
	tmp = polrec(1.0, lat2/cf); rxy=tmp[1]; z2=tmp[2]
	tmp = polrec(rxy, long2/cf); x2=tmp[1]; y2=tmp[2]

	#--- Compute vector dot product for both points. ---
	cs = x1*x2 + y1*y2 + z1*z2

	#--- Compute the vector cross product for both points. ---
	xc = y1*z2 - z1*y2
	yc = z1*x2 - x1*z2
	zc = x1*y2 - y1*x2
	sn = sqrt(xc*xc + yc*yc + zc*zc)

	#--- Convert to polar.  ------
	recpol = function(x,y) c(sqrt(x*x+y*y), atan2(y,x))
	tmp = recpol(cs, sn); r=tmp[1]; a=tmp[2]
	return(cf*a)
 }

