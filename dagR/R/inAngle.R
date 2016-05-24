inAngle <-
function(a,b)
{ # internally used by smoothArc();
  rv<-b-a;
  if(abs(rv)>pi) rv<-sign(rv)*(-1)*(2*pi-abs(rv));
  return(rv);
}

