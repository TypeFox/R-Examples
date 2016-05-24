Newsboy<-function(m,sd,p,c,s=0){ CR=(p-c)/(p-s);
                                 CV<-sd/m;
                                 z<- qnorm(p=CR, mean = 0, sd = 1);
                                 q=m+z*sd;
                                 G<-(p-s)*sd*dnorm(z)
                                 P<-(p-c)*m-G
                                FR<-1-CV*(dnorm(z)-(1-pnorm(z))*z);
                                 SS=z*sd;
                                 options(digits=2)
                                 f=c(Q=q,SS=SS,ExpC=G,ExpP=P,CV=CV,CR=CR,FR=FR,z=z)
                                 return(f)
                                }   
