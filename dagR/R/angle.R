angle <-
function(A,B)
{ # internally used by smoothArc();
  rv<-atan((B[1]-A[1])/(B[2]-A[2]));
  if(B[2]<A[2]) rv<-rv-sign(rv)*pi;
  if(rv==0 && B[2]<A[2]) rv<-pi;

  rv<-(rv+pi)%%(2*pi);    # this changes the range from (-pi, pi);
                          # to (0, 2pi);

  return(rv);
}

