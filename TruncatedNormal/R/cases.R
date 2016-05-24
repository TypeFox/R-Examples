cases <-
  function(p,l,u){
    x=rep(NaN,length(l))
    a=35; # treshold for switching between erfcinv and newton method
    # case 1: a<l<u
    I=l>a;
    if (any(I)){
      tl=l[I]; tu=u[I]; tp=p[I]; x[I]=normq(tp,tl,tu);
    }
    # case 2: l<u<-a
    J=(u<(-a));
    if (any(J)){
      tl=-u[J]; tu=-l[J]; tp=p[J]; x[J]=-normq(1-tp,tl,tu);
    }
    # case 3: otherwise use erfcinv
    I=!(I|J);
    if  (any(I)){
      tl=l[I]; tu=u[I]; tp=p[I]; x[I]=Phinv(tp,tl,tu);
    }
    return(x)
  }
