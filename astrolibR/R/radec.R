radec = function(ra,dec, hours=F) {
  
  if(hours){
    ra = ra %% 24.
    ra = ra + 24*(ra<0)
    ihr = as.integer(ra)
    xmin = abs(ra*60. - ihr*60.)
  } else {
    ra = ra %% 360.          #Make sure between 0 and 24 hours
    ra = ra + 360*(ra<0)
    ihr = as.integer(ra/15.)
    xmin =abs(ra*4.0-ihr*60.0)
  }
  imin = as.integer(xmin)
  xsec = (xmin-imin)*60.0
  ideg = as.integer(dec)
  xmn = abs(dec-ideg)*60.0
  imn = as.integer(xmn)
  xsc = (xmn-imn)*60.0
  zero_deg = ( ideg==0 ) & (dec<0)
  imn = imn - 2*imn*as.integer( zero_deg*(imn!=0) )
  xsc = xsc - 2*xsc*zero_deg*(imn==0)
  return(list(ihr=ihr,imin=imin,xsec=xsec,ideg=ideg,imn=imn,xsc=xsc))
}
