`lowplane` <-
function( az, dip, col = par("col"), UP=FALSE, PLOT=TRUE)
  {
    
    if(missing(col)) { col = par("col")  }
    if(missing(PLOT)) { PLOT = TRUE }
#      az = strike of the plane (NOT down dip azimuth)
#      dip = dip from horizontal
#  lam is the azimuth  
#  phi is the azimuth
    
#    A = MOD(az, 360)
     DEG2RAD = pi/180;
#    lam = A$rem
    beta = az * DEG2RAD;
    
   #  print(paste(" ","LOWPLANE: ", az, dip, UP))
                  #   given the dip and strike of a plan, plot it
    if(missing(UP)) { UP=FALSE }
    
                          #  if(!UP) { dip = 90-dip; lam=lam+180 }
    if(UP==TRUE) { beta = beta+pi }
    
  #   beta = (90+az)* DEG2RAD;
     
    pi180 = pi / 180;
    phi = pi180*seq(-90,90, by=1);
    rz = matrix(ncol=2, nrow=2)
#    rtt = -lam
                                        # rtt = 90-lam
    rz[1,1] = cos(beta);
    rz[1,2] = sin(beta);
    rz[2,1] = -rz[1,2];
    rz[2,2] = rz[1,1];
  co = cos(beta);
  si = sin(beta);
 
    if(dip != 0)
      {
        lambda =   (90-dip) *  DEG2RAD;
        alpha = acos(cos(phi) * cos(lambda));
        tq = sqrt(2)*sin(alpha/2.0);
        sint = sin(phi) / sin(alpha);
        temps = rep(1,length(sint)) -  (sint  * sint)
        temps[is.nan(temps)] = 0
        temps[temps<0] = 0
        x = tq * sqrt(temps  );
        y = tq * sint   ; 
        # prj = cbind(x,y) %*% rz
        
       # x1 = prj[,1]
       # y1 = prj[,2]
        x1 = co * x + si * y;
        y1 =  -si * x + co * y;  
      }
    else
      {

         x = cos(phi)
        y = sin(phi)
         x1 = co * x + si * y;
         y1 =  -si * x + co * y;  
         x = c(0,0)
        y = c(-1,1)
        prj = cbind(x,y) %*% rz 
        x1 = prj[,1]
        y1 = prj[,2]
      }        
    if(PLOT) lines(x1,y1, lwd=2, col=col)
    return(list(x=x1, y=y1))
  }

