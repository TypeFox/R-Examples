euler= function(ai,bi,select, fk4=F,  radian=F) {

  twopi   =   2*pi
  fourpi  =   4*pi
  rad_to_deg = 180/pi
  if(fk4) { 
    equinox = '(b1950)' 
    psi   = c( 0.57595865315, 4.9261918136,  
      0.00000000000, 0.0000000000,    
      0.11129056012, 4.7005372834)
    
    stheta =c(0.88781538514,-0.88781538514, 
      0.39788119938,-0.39788119938, 
      0.86766174755,-0.86766174755)
    
    ctheta =c(0.46019978478, 0.46019978478, 
      0.91743694670, 0.91743694670, 
      0.49715499774, 0.49715499774)
    
    phi  = c(4.9261918136,  0.57595865315, 
      0.0000000000, 0.00000000000, 
      4.7005372834, 0.11129056012)
    
  }
  else { 
    equinox = '(j2000)'
    psi   = c( 0.57477043300, 4.9368292465,  
            0.00000000000, 0.0000000000,    
            0.11142137093, 4.71279419371)
    stheta =c( 0.88998808748,-0.88998808748, 
      0.39777715593,-0.39777715593, 
      0.86766622025,-0.86766622025)
    
    ctheta =c( 0.45598377618, 0.45598377618, 
      0.91748206207, 0.91748206207, 
      0.49714719172, 0.49714719172)
    
    phi  = c( 4.9368292465,  0.57477043300, 
      0.0000000000, 0.00000000000, 
      4.71279419371, 0.11142137093)
    
  }


  i  = select #- 1     (not needed in R) # idl offset

     if(radian) { 
       ao = ai - phi[i]
       bo = bi
     }
     else {      
       ao  = ai/rad_to_deg - phi[i]
       bo = bi/rad_to_deg
     }	  
   
  sb = sin(bo)
  cb = cos(bo)
  cbsa = cb * sin(ao)
  bo  = -stheta[i] * cbsa + ctheta[i] * sb
  tmp = bo
  tmp[tmp>1]=1
  bo    = asin(tmp) 
  if(!radian) bo = bo*rad_to_deg

  ao =  atan2( ctheta[i] * cbsa + stheta[i] * sb, cb * cos(ao) )
  ao = ( (ao+psi[i]+fourpi) %% twopi) 
  if(!radian) ao = ao*rad_to_deg
  browser()
  return(list(ao=ao,bo=bo))
}
