#source('co_nutate.R')
#source('ct2lst.R')
#source('altaz2hadec.R')
#source('precess.R')
#source('co_refract.R')

hor2eq = function( alt, az, jd,  lat=43.0783, lon= -89.865, ws=FALSE,
           b1950 = FALSE, precess_=TRUE, nutate_=TRUE,
           refract_ = TRUE, aberration_ = TRUE, altitude=0) {


    d2r = pi/180.
    h2e = 15

    if (refract_) alt = co_refract(alt, altitude=altitude)

    if(ws) az = az - 180


    tmp = co_nutate(jd, 45.,45.)
    dra1 = tmp$d_ra
    ddec1 = tmp$d_dec
    eps=tmp$eps
    d_psi= tmp$d_psi
        

    lmst = ct2lst(lon, 0, jd)
    lmst = lmst*h2e # convert lmst to degrees (btw, this is the ra of the zenith)
    last = lmst + d_psi *cos(eps)/3600. # add correction in degrees
     

    tmp = altaz2hadec(alt, az, lat)
    ha=tmp$ha
    dec=tmp$dec
    ra = (last - ha + 360.) %% 360.

    tmp = co_nutate( jd, ra, dec)
    dra1 = tmp$d_ra
    ddec1 = tmp$d_dec
    eps=tmp$eps
    d_psi= tmp$d_psi


    tmp = co_aberration( jd, ra, dec, eps)
    dra2 = tmp$d_ra
    ddec2 = tmp$d_dec
    eps=tmp$eps

    ra = ra - (dra1*nutate_ + dra2*aberration_)/3600.
    dec = dec - (ddec1*nutate_ + ddec2*aberration_)/3600.
    j_now = (jd - 2451545.)/365.25 + 2000.0 # compute current equinox
    njd = length(j_now)
    npos = length(ra)
    if ((njd==1) && (npos>1)) j_now = rep(j_now, npos)

    if (precess_) {
        if(b1950){
            for(i in 1:npos) {
                tmp = precess(ra[i], dec[i], j_now[i], 1950.0, fk4=TRUE)
                ra[i] = tmp$ra
                dec[i] = tmp$dec
            }
        } 
	else {
            for(i in 1:npos) {
                tmp = precess(ra[i], dec[i], j_now[i], 2000.0)
                ra[i] = tmp$ra
                dec[i] = tmp$dec
            }
        }
    }

    return(list(ra=ra,dec=dec,ha=ha))
}

