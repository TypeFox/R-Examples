confband<-function(sobj,conf.level=0.95){

km<-sobj

k<-sum(sobj$n.event>0)

bands<-matrix(0,k+1,2)

g<-function(x) cover(x,km) - conf.level

lower<-upper<-iv(km)

cover2<-cover1<-g(lower)

if (cover1<=0) {while (cover2<=0){upper<-upper+0.1
                                  cover2<-g(upper)}}
               else {while (cover1>0) {lower<-lower-0.05
                                       cover1<-g(lower)}}

critval<-zbrent(g,c(lower,upper),1.0e-80)

bands<-exact(km,critval)

cat("The critical value required is ",critval,"\n")

bands}									   
