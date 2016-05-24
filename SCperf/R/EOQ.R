EOQ=function(d,k,h,b=0){ if (b == 0){Q<-(2*d*k/h)^0.5;
                                  T<-Q/d;
                                  TVC<-(2*d*k*h)^0.5;
                                  f<-c(Q=Q,T=T,TVC=TVC)
                                                     
                                  }
                        
                         else{ Q<-(2*d*k*(h+b)/(h*b))^0.5;
                               T<-Q/d; 
                               s<-Q*(h/(h+b));
                              TVC<-(2*k*d*h*b/(h+b))^0.5;                             
                               f<-c(Q=Q,T=T,S=s,TVC=TVC)
                              }
                          options(digits=2,scipen=3)
                         return(f)
                        }
