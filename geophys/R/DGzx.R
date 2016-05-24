DGzx<-function(xs, zs,  xv, zv, den)
  {
  ###  Compute the gravity anomaly following Won and Bevis
    twopi = 2*pi
    con=13.3464E-03 

    nvert = length(xv)

    ####  make sure the polygon is closed

    if(xv[1]!=xv[nvert] & zv[1]!=zv[nvert])
      {
        xv = c(xv, xv[1])
        zv = c(zv, zv[1])
        
      }
    nvert = length(xv)
###  create output vector
    gravz = rep(NA, length(xs))
    gravx = rep(NA, length(xs))
    ###  loop over the stations
    for(i in 1:length(xs))
      {

        xst = xs[i]
        zst = zs[i]
        
        x1 = xv[1:(nvert-1)]-xst;
        z1 = zv[1:(nvert-1)]-zst;
        x2 = xv[2:(nvert)]-xst;
        z2 = zv[2:(nvert)]-zst;

        ### calculate the angles
        theta1 = atan2(z1, x1);
        theta2 = atan2(z2, x2);
###  need to get rid of pathology here.

        dsign = sign(z1) != sign(z2)
        if(any(dsign))
          {
        theta1[ dsign & (x1*z2<x2*z1) & z2>=0] = theta1[ dsign & (x1*z2<x2*z1) & z2>=0]+twopi
        theta2[ dsign & (x1*z2>x2*z1) & z1>=0] = theta2[ dsign & (x1*z2>x2*z1) & z1>=0 ]+twopi
      }

        dx = x2-x1;
        dz = z2-z1;

        r1 = sqrt(x1*x1 + z1*z1);
        r2 = sqrt(x2*x2 + z2*z2);

        
        dxz2 = (dx*dx + dz*dz) 
        A = dx*( x1*z2 - x2*z1 )/(dx*dx + dz*dz);

        B = dz/dx;


        ZEE =  A*( (theta1 - theta2)  + B*log(r2/r1))
        EX =  A*(-((theta1 - theta2) )*B + log(r2/r1))

        ZEE[x1*z2 == x2*z1] = 0
        EX[x1*z2 == x2*z1] = 0

        ZEE[ (x1==0 & z1==0) |  (x2==0 & z2==0)   ] = 0
        EX[ (x1==0 & z1==0) |  (x2==0 & z2==0)   ] = 0

        ZEE[x1==x2] = x1[x1==x2] * log(r2[x1==x2]/r1[x1==x2])
        
        EX[x1==x2] = -1*x1[x1==x2] *(theta1[x1==x2] - theta2[x1==x2])
        
        Z = sum( ZEE );
        X = sum( EX );

        gravz[i] = con*den*Z
        gravx[i] = con*den*X
        
      }

    invisible(list(Gz=gravz, Gx=gravx))
    
}
