EPQ=function(d,p,k,h){nn<-2*d*k;
                      dd<-h*(1-d/p) ;
                      Q<-(nn/dd)^0.5;
                      t<-Q/p
                      T<-Q/d;
                      I<-(1-d/p)*Q
                      TC<-k*(d/Q)+h*(Q/2)*(1-d/p);                             
                      f<-c(q=Q,t=t,T=T,I=I,TC=TC)
                         
                      options(digits=6)  
                      return(f)
                      }


