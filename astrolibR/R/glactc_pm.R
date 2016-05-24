#source('glactc.R')


glactc_pm = function(ra,dec,mu_ra,mu_dec,year,gl,gb,mu_gl,mu_gb,j,
	degree=FALSE, fk4 = FALSE, supergalactic = FALSE, mustar=FALSE) {

  radeg = 180.0/pi
  if(supergalactic ){
    rapol = 283.18940711/15.0
    decpol = 15.64407736
    dlon =  26.73153707
  }
  else {
    rapol = 12.0e0 + 49.0e0/60.0e0
    decpol = 27.4e0
    dlon = 123.0e0
  }
  
  sdp = sin(decpol/radeg)
  cdp = sqrt(1.0e0-sdp*sdp)
  radhrs=radeg/15.0e0
 
if(j==1) {
       if(!degree)
         ras = ra*15.0
       else
         ras =ra

       decs = dec
       if(!fk4){
         if(year!=2000 ) {
           tmp = precess(ras,decs,year,2000)
           ras = tmp$ra
           decs = tmp$dec
         }
         tmp = bprecess(ras,decs)
         ras = tmp$ra_1950
         decs = tmp$dec_1950
       }
       else if(year!=1950 ) {
         tmp = precess(ras,decs,year,1950,fk4=T)
         ras = tmp$ra
         decs = tmp$dec
       }
       
       raindeg = ras
       ras = ras/radeg - rapol/radhrs
       sdec = sin(decs/radeg)
       cdec = sqrt(1.0e0-sdec*sdec)
       sgb = sdec*sdp + cdec*cdp*cos(ras)
       gb = radeg * asin(sgb)
       cgb = sqrt(1.0e0-sgb*sgb)
       sine = cdec * sin(ras) / cgb
       cose = (sdec-sdp*sgb) / (cdp*cgb)
       gl = dlon - radeg*atan2(sine,cose)
       ltzero=(gl<0.0)
       gl[ltzero]=gl[ltzero]+360.0e0
       if(!mustar)mu_ra = mu_ra*cdec
       mu_gb =
         mu_dec*(cdec*sdp-sdec*cdp*cos(ras))/cgb - mu_ra*cdp*sin(ras)/cgb
       mu_gl = sqrt(mu_dec^2 + mu_ra^2 - mu_gb^2)
       tmp = glactc(raindeg, decs, year, gl0, gb0,1,
         degree=TRUE,supergalactic=supergalactic)

       gl0 = tmp$gl
       gb0 = tmp$gb
       radelta = 1e-2*mu_ra/abs(mu_ra)
       decdelta = decs + 1e-2*mu_dec/abs(mu_ra)
       tmp = glactc(raindeg+radelta, decs+decdelta, year, gl2, gb2,
              1,degree=TRUE,supergalactic=supergalactic)
       gl2 = tmp$gl
       gb2 = tmp$gb

       if (gl2<gl0) mu_gl = -mu_gl
       if(!mustar)mu_gl = mu_gl/cgb
        return(list(gl=gl,gb=gb,mu_gl= mu_gl, mu_gb = mu_gb))
     }
    else if(j==2)  {
        sgb = sin(gb/radeg)
        cgb = sqrt(1.0e0-sgb*sgb)
        sdec = sgb*sdp + cgb*cdp*cos((dlon-gl)/radeg)
        dec = radeg * asin(sdec)
        cdec = sqrt(1.0e0-sdec*sdec)
        sinf = cgb * sin((dlon-gl)/radeg) / cdec
        cosf = (sgb-sdp*sdec) / (cdp*cdec)
        ra = rapol + radhrs*atan2(sinf,cosf)
        ra = ra*15.0
        if(!mustar) mu_gl = mu_gl*cgb
        mu_dec =
          mu_gb*(cgb*sdp-sgb*cdp*cos((dlon-gl)/radeg))/cdec+
            mu_gl*cdp*sin((dlon-gl)/radeg)/cdec
        mu_ra = sqrt(mu_gl^2 + mu_gb^2 - mu_dec^2)
        mu_gl_delta = 1e-2*mu_gl/abs(mu_gl)
        mu_gb_delta = 1e-2*mu_gb/abs(mu_gl)
        tmp = glactc(ra2, dec2, 1950, gl+mu_gl_delta, gb+mu_gb_delta,
               2,degree=TRUE,supergalactic=supergalactic)
        ra2 = tmp$ra
        dec2 = tmp$dec
        if (ra2<ra) mu_ra = -mu_ra
  if(!mustar) mu_ra = mu_ra/cdec
        if(!fk4) {
          ras = ra
          decs = dec
          tmp = jprecess(ras,decs)
          ra = tmp$ra_2000
          dec = tmp$dec_2000
          
          if(year!=2000 ) {
            tmp = precess(ra,dec,2000,year)
            ra = tmp$ra
            dec = tmp$dec
          }
        }
        else if(year!=1950) {
          tmp = precess(ra,dec,1950,year,fk4=TRUE)
          ra = tmp$ra
          dec = tmp$dec
        }

        gt36 = (ra>360.0)
       ra[gt36] = ra[gt36] - 360.0e0
       if(!degree)     ra = ra / 15.0e0
        return(list(ra=ra, dec=dec, mu_ra = mu_ra, mu_dec = mu_dec))
      }
}
