#source('ct2lst.R')
#source('hadec2altaz.R')
#source('co_nutate.R')
#source('co_refract.R')
#source('co_aberration.R')
#source('precess.R')

eq2hor = function(
  ra, dec, jd,
  lat= 43.0783, lon= -89.865, ws=F, obsname,
  b1950, precess_=TRUE, nutate_=TRUE,
  refract_=TRUE, aberration_=TRUE,
  altitude=0, ...) {

                                        # 43.0783 is the declination of the zenith

  d2r = pi/180.
  h2r = pi/12.
  h2e = 15

  if (!missing(b1950) )
    s_now='   (j1950)'
  else
    s_now='   (j2000)'

  j_now = (jd - 2451545.)/365.25 + 2000.0 # compute current equinox
  if (precess_ ){
      if(!missing(b1950) ){
          for(i in 1:length(jd)) {
              tmp=precess(ra[i], dec[i], 1950.0, j_now[i], fk4=TRUE)
              ra[i] = tmp$ra
              dec[i] = tmp$dec
          }
      }
      else {
          for(i in 1:length(jd)) {
              tmp=precess(ra[i], dec[i], 2000.0, j_now[i])
              ra[i] = tmp$ra
              dec[i] = tmp$dec
          }
      }

      tmp = co_nutate(jd, ra, dec)
      dra1 = tmp$d_ra
      ddec1 = tmp$d_dec
      eps=tmp$eps
      d_psi=tmp$d_psi

      tmp = co_aberration(jd, ra, dec,eps)
      dra2 = tmp$d_ra
      ddec2 = tmp$d_dec
      eps=tmp$eps

      ra = ra + (dra1*nutate_ + dra2*aberration_)/3600.
      dec = dec + (ddec1*nutate_ + ddec2*aberration_)/3600.

      lmst = ct2lst(lon, 0, jd)  # get lst (in hours) - note:this is independent of
                                        #time zone  since giving jd
      lmst = lmst*h2e # convert lmst to degrees (btw, this is the ra of the zenith)
      last = lmst + d_psi *cos(eps)/3600. # add correction in degrees

      ha = last - ra
      w = (ha<0)
      ha[w] = ha[w] + 360.
      ha = ha %% 360

      tmp = hadec2altaz(ha, dec, lat, ws=ws)
      alt = tmp$alt
      az = tmp$az

      if (refract_ )
          alt = co_refract(alt, altitude=altitude, ..., to_observed=TRUE)
  }
  return(list(alt=alt, az=az, ha=ha))
}
