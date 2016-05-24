`polyfoc` <-
function(strike1, dip1, strike2, dip2, PLOT=FALSE, UP=TRUE)
  {
    if(missing(UP)) {  UP=TRUE }
    if(missing(PLOT)) {  PLOT=TRUE }
    RAD2DEG = 180/pi;
    num = 40;
    F1 = faultplane(strike1, dip1,  PLOT=FALSE, UP=UP)
    F2 = faultplane(strike2, dip2,   PLOT=FALSE,  UP=UP)

   
    k = length(F1$x)
    k2 = length(F2$x)

    ################ here the angle depends on the direction taken by the plane plot
    ####   must be careful

    
    alpha1 = atan2(F1$y[1], F1$x[1]);
    alpha2 = atan2(F1$y[k], F1$x[k]);
    
    beta1 = atan2(F2$y[1], F2$x[1]);
    beta2 = atan2(F2$y[k2], F2$x[k2]);


    a1 = RAD2DEG * alpha1;
    a2 = RAD2DEG * alpha2;
    b1 = RAD2DEG * beta1;
    b2 = RAD2DEG * beta2;

    tx1= rep(0,num);
    ty1  = rep(0,num);
    
 
    dang = bang(F1$x[k], F1$y[k], F2$x[1], F2$y[1])/num

   
    ang = alpha2+dang*seq(1,num);
    tx1  = cos(ang);
    ty1  = sin(ang);
     
    
    tx2= rep(0,num);
    ty2  = rep(0,num);
    
    dang = bang(F2$x[k2], F2$y[k2], F1$x[1], F1$y[1])/num
    
    ang = beta2+dang*seq(1,num) ;
    
    ang = ang+dang;
    tx2 = cos(ang);
    ty2 = sin(ang);
    
    
    Px = c(F1$x, tx1, F2$x, tx2);
    Py = c(F1$y, ty1, F2$y, ty2);

    if(PLOT==TRUE)
      {
        lines(Px, Py)
      }
    
    return(list(Px=Px, Py=Py))
  }

