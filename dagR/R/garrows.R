garrows <-
function(x0,y0,x1,y1, xgap, ygap, len=0.1)
{ # draws the arrows in the DAG...;
  if((x0<x1)&&(abs(atan((y1-y0)/(x1-x0))*180/pi)<60))
  { xx0<-x0+xgap;
    xx1<-x1-xgap;
  } else if ((x0>x1)&&(abs(atan((y1-y0)/(x1-x0))*180/pi)<60))
  { xx0<-x0-xgap;
    xx1<-x1+xgap;
  } else
  { xx0<-x0;
    xx1<-x1;
  } 
  if((y0<y1)&&(abs(atan((y1-y0)/(x1-x0))*180/pi)>30))
  { yy0<-y0+ygap;
    yy1<-y1-ygap;
  } else if ((y0>y1)&&(abs(atan((y1-y0)/(x1-x0))*180/pi)>30))
  { yy0<-y0-ygap;
    yy1<-y1+ygap;
  } else
  { yy0<-y0;
    yy1<-y1;
  }
  arrows(xx0, yy0, xx1, yy1, length=len);
}

