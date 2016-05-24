ShadowCLVD <-
function(m, PLOT=TRUE, col=rgb(1, .75, .75))
{### start ShadowCLVD

### -----------------------------------------
###  lower hemisphere equal-area projection
### -----------------------------------------
###   originally developed for matlab by Vaclav Vavrycuk
### see: http://www.ig.cas.cz/en/research-&-teaching/software-download/
  ### translated to R by Jonathan M. Lees

  if(missing(col)) { col =rgb(1, .75, .75)  }
  if(missing(PLOT)) { PLOT=TRUE }
  
  
  x.min = -1
  x.max = 1
  dx = .025

  y.min = -1
  y.max = 1
  dy = .025

  
  EX = seq(from=x.min, by=dx, to=x.max)
  WHY = seq(from=y.min, by=dy, to=y.max)
  ndim = length(EX)
  
  M = RPMG::meshgrid(EX ,WHY  )

  
  r = sqrt(M$x^2+M$y^2)

  sin.fi = M$x/r
  cos.fi = M$y/r
  
  sin.fi[r<=1.e-5]  = 0
  cos.fi[r<=1.e-5] = 0


  theta = matrix(NA, ncol=ncol(r), nrow=nrow(r) )

  theta[r<1] = asin(sqrt((M$x[r<1]^2+M$y[r<1]^2)/2))*360/pi

  n = array(data = NA, dim = c(ncol(r),nrow(r), 3 )  , dimnames = NULL)


  n[,,2]  = sin(theta*pi/180)*sin.fi   ###  n[2] directed to the East
  n[,,1]  = sin(theta*pi/180)*cos.fi   ###  n[1] directed to the North
  n[,,3]  = sqrt(1-n[,,1]^2-n[,,2]^2)

  u.radiation.z = matrix(0, ncol=ncol(r), nrow=nrow(r) )

  for (  i  in  1:3 )
    {
      for (  j  in  1:3 )
        {
          u.radiation.z =  u.radiation.z + n[,,3]*n[,,i]*n[,,j]*m[i,j]
        }
    }
  
  sign.u.radiation.z = sign(u.radiation.z)


  
  u.radiation.z[ sign.u.radiation.z<=0]  = NA

z1 = t(u.radiation.z)

  
z2 = z1
  

  if(PLOT)  image(x=EX, y=WHY, z=z2, col=col  , add=TRUE)


  
  invisible( list(x=EX, y=WHY, z=z2) )

}
