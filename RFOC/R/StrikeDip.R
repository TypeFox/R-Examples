StrikeDip<-function(x=x, y=y, MEC, focsiz, addDIP=TRUE, ...)
  {
###  Author: Jonathan Lees, Univerity of North Carolina, Feb 2010
                                        #  input:
###    x,y = location on the plot
###  MEC is theoutput of SDRfoc  program
###    focsiz  size of the focal sphere
    
#### using a strike and dip of a focal mechanism
#### plot a segment along the strike and
####a perpendicular segment representing the dip
###    see geotouch: this is the ralph archuletta style

###   Plot Style: planes(only)
###     planes/PT
###     PT (only)
###     add fault plane
###     fault plane only

    ##  points(x, y, col='blue', pch=3)

    if(missing(addDIP)) { addDIP=TRUE }
    
    strk1 = MEC$az1
    dip1 = MEC$dip1
    
    tq = 1.0;
    xi =  0.0174532 * strk1 - 1.570796327;
    ang = strk1*pi/180
    dip  = dip1*pi/180
    
    p1x = x + focsiz*sin(ang)
    p1y = y + focsiz*cos(ang)
    p2x = x - focsiz*sin(ang)
    p2y = y - focsiz*cos(ang)
    
### points(p1x, p1y, col='blue', pch=3)
###  points(p2x, p2y, col='blue', pch=3)

    segments(p1x, p1y, p2x, p2y, ...)
    
    
####	 XDrawLine(dpy, xwin, gc, kx, ky, kx2, ky2);
    
####	 /*   Down  Dip line    */

    if(addDIP)
      {
        F = focpoint(MEC$F$az, MEC$F$dip, col = 'black' , lab = "F", UP = MEC$UP, PLOT = FALSE)


        p3x = x + focsiz*F$x
        p3y = y + focsiz*F$y

        segments(x, y, p3x, p3y, ...)

      }


  }
