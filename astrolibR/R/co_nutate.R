#source('nutate.R')

co_nutate = function(jd, ra, dec) {

  d2r = pi/180
  d2as = pi/(180*3600)
  T = (jd -2451545.0)/36525.0 # Julian centuries from J2000 of jd.
  tmp=nutate(jd)
  d_psi=tmp$nut_long
  d_eps = tmp$nut_obliq
  
  eps0 = 23.4392911*3600. - 46.8150*T - 0.00059*T^2 + 0.001813*T^3
  eps = (eps0 + d_eps)/3600.*d2r # true obliquity of the ecliptic in radians
  ce = cos(eps)
  se = sin(eps)
  x = cos(ra*d2r) * cos(dec*d2r)
  y = sin(ra*d2r) * cos(dec*d2r)
  z = sin(dec*d2r)
  x2 = x - (y*ce + z*se)*d_psi * d2as
  y2 = y + (x*ce*d_psi - z*d_eps) * d2as
  z2 = z + (x*se*d_psi + y*d_eps) * d2as
  r = sqrt(x2^2 + y2^2 + z2^2)
  xyproj = sqrt(x2^2 + y2^2)
  ra2 = x2 * 0.
  dec2= x2 * 0.
  w1 = ( (xyproj==0) & (z!=0) )
  w2 = (xyproj!=0)
                                        # places where xyproj=0 (point at NCP or SCP)
  dec2[w1] = asin(z2[w1]/r[w1])
  ra2[w1] = 0.
  
                                        # places other than NCP or SCP
  ra2[w2] = atan2(y2[w2],x2[w2])
  dec2[w2] = asin(z2[w2]/r[w2])
                                        # convert to DEGREES
  ra2 = ra2 /d2r
  dec2 = dec2 /d2r
  w = which(ra2<0.)
  ra2[w] = ra2[w] + 360.
  d_ra = (ra2 - ra) * 3600.
  d_dec = (dec2 - dec) * 3600.

  return(list(d_ra=d_ra, d_dec=d_dec,
              eps=eps, d_psi=d_psi, d_eps=d_eps))  
}
