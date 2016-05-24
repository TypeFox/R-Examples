adstring = function(ra_dec,dec,precision, truncate = FALSE) {
  
  if(missing(dec)) {
    stopifnot(length(ra_dec)==2)
    ra = ra_dec[1] %% 360
    dec = ra_dec[2]
  }
  else { 	
   	stopifnot(length(ra_dec)==length(dec))
  	ra = ra_dec %% 360.
  }

  
  badrange = (dec< -90) | (dec > 90)
  if(any(badrange))
    stop('WARNING - Some declination values are out of valid range (-90 < dec <90)')
  tmp = radec(ra, dec)
  ihr = tmp$ihr
  imin = tmp$imin
  xsec = tmp$xsec
  
  ideg = tmp$ideg
  imn = tmp$imn
  xsc = tmp$xsc
  
  if(length(precision)==0) precision = 0
  if(precision<0) precision = 0
  if(precision>4) precision = 4         #No more than 4
                                        #decimal places

                                        #-------------------------------------------

  #Deal with imin, xsec
  if(!truncate) {
    roundsec = c(59.5,59.95,59.995,59.9995,59.99995,59.999995)
    carry = which(xsec>roundsec[precision+2])
    imin[carry] = imin[carry] + 1
    xsec[carry] = 0.0
    
    mcarry = which(imin[carry]==60)
    ic = carry[mcarry]
    ihr[ic] = (ihr[ic] + 1) %% 24
    imin[ic] = 0
  }
  else {
    xsec = trunc(xsec,precision)
  }
  
  secfmt = sprintf('%%%d.%df',precision+4,precision+1)

  fullfmt = paste('%2.0f%3.0f ',secfmt,sep='')
  result = sprintf(fullfmt,ihr,imin,xsec)
  
  if(length(precision)==0) precision = 1

  #Deal with imn, xsc
  imn = abs(imn)
  xsc = abs(xsc)
  if(( precision==0 ) ){ 
    secfmt = '%2.0f' 
    if(!truncate){
      xsc = round(xsc)
      #Carry sec into min
      carry = (xsc==60)
      if(any(carry)){                 #Updated April 2002
        xsc[carry] = 0
        imn[carry] = imn[carry] + 1
      }
    }
  }
  else {
      secfmt = sprintf('%%%d.%df',precision+3,precision)
    
      if(!truncate){
        ixsc = as.integer(xsc + 0.5/10^precision)
                                        #Carry sec into min
        carry = ixsc>=60
        xsc[carry] = 0
        imn[carry] = imn[carry] + 1
      }
      else {
        xsc = trunc(xsc,precision)
      }
    }

  #Carry min into deg
  
  pos = dec>=0 
  carry = (imn==60); 
  if(any(carry)){
    ideg[carry] = ideg[carry] - 1 + 2*pos[carry]
    imn[carry] = 0
  }
  deg = sprintf('%+3.2d',ideg)

  #Adjust deg for sign
  zero = (ideg==0)
  if(any(zero)){
    negzero = (dec[zero]<0);
    if(any(negzero)) deg[zero[negzero]] = '-00' 
  }

  paste(result,' ',deg,sprintf('%3.2d',imn),' ',sprintf(secfmt,xsc),sep='')
}
