altaz2hadec=function( alt, az, lat) {

 d2r = pi/180
 alt_r  = alt*d2r
 az_r = az*d2r
 lat_r = lat*d2r
 ha = atan2( -sin(az_r)*cos(alt_r), 
           -cos(az_r)*sin(lat_r)*cos(alt_r)+sin(alt_r)*cos(lat_r))
 ha = ha / d2r
 w = which(ha<0)
 ha[w] = ha[w] + 360.
 ha = ha %% 360
 sindec = sin(lat_r)*sin(alt_r) + cos(lat_r)*cos(alt_r)*cos(az_r)
 dec = asin(sindec)/d2r  # convert dec to degrees
 return(list(ha=ha, dec=dec))
}
