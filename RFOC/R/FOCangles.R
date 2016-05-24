FOCangles <-
function( m)
{
### *   function calculates a fault normal and slip from the moment tensor
### *                                                        

###   originally developed for matlab by Vaclav Vavrycuk
### see: http://www.ig.cas.cz/en/research-&-teaching/software-download/
### translated to R by Jonathan M. Lees
  
### ---------------------------------------------------------
###  eigenvalues and eigenvectors of a moment tensor
### ---------------------------------------------------------
  Vnorm <-function(X){ return( sqrt(sum(X^2)) ) }

  JEIG = eigen(m)
  nv = length(JEIG$values)
#########################  matlab does the reverse order
  vector = JEIG$vectors[,rev(1:nv)  ] 
  diag.m= rev(JEIG$values)

  values.total = diag( diag.m )
  

### -----------------------------------------------------
###  isotropic part of the moment tensor
### ----------------------------------------------------
  volumetric = (1/3)*sum(diag.m)

  m.volumetric = volumetric*diag(1,3)

### ----------------------------------------------------
###  deviatoric part of the moment tensor
### ---------------------------------------------------
  m.deviatoric = diag(diag.m)  - m.volumetric 
### ---------------------------------------------------
###  eigenvalues of the deviatoric part of the moment tensor
### --------------------------------------------------

  Edev = rev( eigen(m.deviatoric, only.values =TRUE)$values )

  
  value.dev = Edev
  
### --------------------------------------------------
###  calculation of a fault normal and slip from the moment tensor
### -------------------------------------------------
  j=order(value.dev)
  values=value.dev[j]
  
  n1 = vector[,j[3]]+vector[,j[1]]
  
  n1 = n1/Vnorm(n1)
  
  if (n1[3]>0) { n1 = -n1    }
###  vertical component is always negative!
  
  u1 = vector[,j[3]]-vector[,j[1]]
  
  u1 = u1/Vnorm(u1)

  if (( n1%*% m %*% u1 + u1 %*% m %*% n1 ) < 0) { u1 = -u1    }
  

  n2 = u1
  u2 = n1
  
  if (n2[3]>0){ n2 = -n2;  u2 = -u2 }
                                        #   ###  vertical component is always negative!
### -------------------------------------------------------
###  1st solution
### ----------------------------------------------------------
  dip    = acos(-n1[3])*180/pi
  strike = asin(-n1[1]/sqrt(n1[1]^2+n1[2]^2))*180/pi

###  determination of a quadrant
  if (n1[2]<0) { strike=180-strike    } #

  rake = asin(-u1[3]/sin(dip*pi/180))*180/pi

###  determination of a quadrant
  cos.rake = u1[1]*cos(strike*pi/180)+u1[2]*sin(strike*pi/180)
  if (cos.rake<0) { rake=180-rake    } #

  if(strike < 0   ) { strike = strike+360    } 
  if(rake  < (-180) ){  rake   = rake  +360    } 
  if(rake  > 180) { rake   = rake  -360    }  #  rake is in the interval -180<rake<180
  
  angles1 = c(strike, dip, rake)

### ---------------------------------------------
###  2nd solution
### -----------------------------------------------
  dip    = acos(-n2[3])*180/pi
  strike = asin(-n2[1]/sqrt(n2[1]^2+n2[2]^2))*180/pi

###  determination of a quadrant
  if (n2[2]<0) {  strike=180-strike    } 

  rake = asin(-u2[3]/sin(dip*pi/180))*180/pi

###  determination of a quadrant
  cos.rake = u2[1]*cos(strike*pi/180)+u2[2]*sin(strike*pi/180)
  if (cos.rake < 0) { rake=180-rake    } 

  if (strike < 0   ){  strike = strike+360    } 
  if (rake  < (-180) ){ rake   = rake  +360    } 
  if (rake  > 180){ rake   = rake  -360    }   
  
  angles2 = c(strike, dip, rake)
### -------------------------------------
  a = c(angles1, angles2)

  return(a)

}
