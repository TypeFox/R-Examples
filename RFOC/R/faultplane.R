`faultplane` <-
function(az,  dip,  col = par("col"), PLOT=TRUE, UP=FALSE) 
  {    
#      az = strike of the plane (NOT down dip azimuth)
#      dip = dip from horizontal
    #   given the dip and strike of a plan, plot it
    DEG2RAD = pi/180;
    beta = az * DEG2RAD;
    
    # print(paste(" ","FaultPLANE: ", az, dip, UP))
                                        
    if(missing(PLOT)) { PLOT=TRUE }
    if(missing(PLOT)) { UP=FALSE }

    
    if(UP==TRUE) { beta = beta+pi }
         
    pi180 = pi / 180;
    phi = pi180*seq(-90,90, by=1);

  co = cos(beta);
  si = sin(beta);
 
    if(dip != 0)
      {
        lambda =   (90-dip) *  DEG2RAD;
        alpha = acos(cos(phi) * cos(lambda));
        tq = sqrt(2)*sin(alpha/2.0);
        
        sint = sin(phi) / sin(alpha);
        sint[is.nan(sint)] = 1
        temps = rep(1,length(sint)) -  (sint  * sint)
        temps[is.nan(temps)] = 0
        temps[temps<0] = 0
        x = tq * sqrt(temps  );
        y = tq * sint   ; 

      }
    else
      {
      ###  x = c(0,0)
      ###  y = c(-1,1)

        x = cos(phi)
        y = sin(phi)
        
       
      }

        x1 =  co   * x +  si * y;
        y1 =  -si  * x +  co * y;
    if(PLOT)
      {
        lines(x1,y1, lwd=2, col=col, xpd=TRUE)
      }
    return(list(x=x1, y=y1))

  }

