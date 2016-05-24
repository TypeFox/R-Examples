#source('precess.R')
#source('precess_xyz.R')

xyz = function(date,equinox) {

  picon = pi/180
  t = (date - 15020.0e0)/36525.0e0         #relative julian century from 1900
  pp = (1.396041 + 0.000308*(t + 0.5))*(t-0.499998)
  el = 279.696678 + 36000.76892*t + 0.000303*t*t - pp
  c = 270.434164 + 480960.*t + 307.883142*t - 0.001133*t*t - pp
  n = 259.183275 - 1800.*t - 134.142008*t + 0.002078*t*t - pp
  g = 358.475833 + 35999.04975*t - 0.00015*t*t
  j = 225.444651 + 2880.0*t + 154.906654*t*t
  v = 212.603219 + 58320.*t + 197.803875*t + 0.001286*t*t
  m = 319.529425 + 19080.*t + 59.8585*t + 0.000181*t*t
  el = el*picon
  g  = g*picon
  j =  j*picon
  c  = c*picon
  v  = v*picon
  n  = n*picon
  m  = m*picon
  x =   0.999860*cos(el)                          
  - 0.025127*cos(g - el)                      
  + 0.008374*cos(g + el)                      
  + 0.000105*cos(g + g + el)                  
  + 0.000063*t*cos(g - el)                    
  + 0.000035*cos(g + g - el)                  
  - 0.000026*sin(g - el - j)                  
  - 0.000021*t*cos(g + el)                    
  + 0.000018*sin(2.*g + el - 2.*v)          
  + 0.000017*cos(c)                           
  - 0.000014*cos(c - 2.*el)                  
  + 0.000012*cos(4.*g + el - 8.*m + 3.*j)  
  - 0.000012*cos(4.*g - el - 8.*m + 3.*j)  
  - 0.000012*cos(g + el - v)                  
  + 0.000011*cos(2.*g + el - 2.*v)          
  + 0.000011*cos(2.*g - el - 2.*j)         
  
  #cat('el=',el,'\n')
  #cat('g=',g,'\n')
  #cat('j=',j,'\n')
  #cat('v=',v,'\n')
  #cat('m=',m,'\n')
  #cat('x=',x,'\n')

  y =   0.917308*sin(el)                             
  + 0.023053*sin(g - el)                         
  + 0.007683*sin(g + el)                         
  + 0.000097*sin(g + g + el)                     
  - 0.000057*t*sin(g - el)                       
  - 0.000032*sin(g + g - el)                     
  - 0.000024*cos(g - el - j)                     
  - 0.000019*t*sin(g + el)                       
  - 0.000017*cos(2.e0*g + el - 2.e0*v)           
  + 0.000016*sin(c)                              
  + 0.000013*sin(c - 2.e0*el )                   
  + 0.000011*sin(4.e0*g + el - 8.e0*m + 3.e0*j)  
  + 0.000011*sin(4.e0*g - el - 8.e0*m + 3.e0*j)  
  - 0.000011*sin(g + el - v)                     
  + 0.000010*sin(2.e0*g + el - 2.e0*v )          
  - 0.000010*sin(2.e0*g - el - 2.e0*j )         
  
  z =   0.397825*sin(el)        
  + 0.009998*sin(g-el)      
  + 0.003332*sin(g+el)      
  + 0.000042*sin(g+g+el)    
  - 0.000025*t*sin(g-el)    
  - 0.000014*sin(g+g-el)    
  - 0.000010*cos(g-el-j)    

  
  if(!missing(equinox)) {
    tmp = precess_xyz( x, y, z, 1950, equinox)
    x = tmp$x
    y = tmp$y
    z = tmp$z
  }
  
  xvel = -0.017200 * sin(el)           
  -0.000288 * sin(g + el)       
  -0.000005 * sin(2.e0*g + el)  
  -0.000004 * sin(c)            
  +0.000003 * sin(c - 2.e0*el)  
  +0.000001 *t * sin(g+el)      
  -0.000001 * sin(2.e0*g-el)           
  
  yvel =  0.015780 * cos(el)            
  +0.000264 * cos(g + el)        
  +0.000005 * cos(2.e0*g + el)   
  +0.000004 * cos(c)             
  +0.000003 * cos(c - 2.e0*el)   
  -0.000001 * t * cos(g + el)    
  zvel = 0.006843 * cos(el)             
  +0.000115 * cos(g  + el)        
  +0.000002 * cos(2.e0*g + el)    
  +0.000002 * cos(c)              
  +0.000001 * cos(c - 2.e0*el)    
  if(!missing(equinox)) {
    tmp = precess_xyz( xvel, yvel, zvel, 1950, equinox)
    xvel = tmp$x
    yvel = tmp$y
    zvel = tmp$z
  }
  return(list(x=x,y=y,z=z,xval=xvel,yval=yvel,zvel=zvel))
}
