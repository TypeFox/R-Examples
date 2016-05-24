sleek <-
function(y)
{
  h<-function(y)
  {
    
  N<-length(y);
  # To calculate runing median length=4 series in z
  z<-NULL;
  z[1]<-y[1];
  z[2]<-median(c(y[1],y[2]));
  z[3]<-median(c(y[2],y[3]));
  
  for(i in 4:(N))
  {
    z[i]<-median(c(y[i-3],y[i-2],y[i-1],y[i]));
  }
  
  z[N-1]<-median(c(y[N-2],y[N-1]));
  z[N]<-median(c(y[N-1],y[N]));
  z[N+1]<-y[N];
  
  # To calculate runing median length=2 series in z1
  z1<-NULL;
  for(i in 1:N)
  {
    z1[i]<-(z[i]+z[i+1])/2;
  }
  
  # To calculate runing median length=5 series in z2
  z2<-NULL;
  z2[1]<-z1[1];
  z2[2]<-median(c(z1[1],z1[2],z1[3]));
  
  for(i in 3:(N-2))
  {
    z2[i]<-median(c(z1[i-2],z1[i-1],z1[i],z1[i+1],z1[i+2]));
  }
  
  z2[N-1]<-median(c(z1[N-2],z1[N-1],z1[N]));
  z2[N]<-z1[N];
  
  # To calculate runing median length=3 series in z3
  z3<-NULL;
  z3[1]<-z2[1];
  
  for(i in 2:(N-1))
  {
    z3[i]<-median(c(z2[i-1],z2[i],z2[i+1]));
  }
  
  z3[N]<-z2[N];
  
  # To calculate Hanning smoother in z4
  z4<-NULL;
  z4[1]<-z3[1];
  
  for(i in 2:(N-1))
  {
    z4[i]<-(z3[i-1]+z3[i]+z3[i+1])/4;
  }
  z4[N]<-z3[N];
  
  # To calculate Endpoint rule
  z4[1]<-median(c(z4[1],z4[2],(3*z4[2]-2*z4[3])));
  z4[N]<-median(c(z4[N],z4[N-1],(3*z4[N-2]-2*z4[N-1]) ));
  
  return(z4)
    }
  # end of 4253H
  
  # To calculate Twicing
  sm<-h(y);
  rf<-(y-sm);
  sm.rf<-h(rf);
  smooth<-(sm.rf+sm);
  
  if(is.ts(y)==1)
  {
    date<-start(y);
    f1<-frequency(y)
    smooth<-ts(smooth, start= date, frequency= f1);
  }
  
  return(smooth)
}
