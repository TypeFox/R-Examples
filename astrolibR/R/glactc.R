#source('precess.R')
#source('bprecess.R')
#source('jprecess.R')

glactc = function(
  ra,dec,year,gl,gb,j, 
  degree= FALSE, fk4= FALSE, supergalactic= FALSE) {

  radeg = 180.0/pi

  if(supergalactic){
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
                                        #
  if(j==1) {
    if(!degree )
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
    return(list(gl=gl,gb=gb))
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
    if(!fk4){
      ras = ra
      decs = dec
      tmp = jprecess(ra,dec)
      ra = tmp$ra_2000
      dec = tmp$dec_2000
      if(year!=2000 ) {
        tmp = precess(ra,dec,2000,year)
        ra = tmp$ra
        dec = tmp$dec
      }
    }
    else if(year!=1950 ) {
        tmp = precess(ra,dec,1950,year,fk4=T)
        ra = tmp$ra
        dec = tmp$dec
    }
  
    
    gt36 =which(ra>360.0)
    ra[gt36] = ra[gt36] - 360.0e0
    if(!degree)     ra = ra / 15.0e0
    return(list(ra=ra,dec=dec))
  }
}
