Fillv=function(nfac,des)
 {
# creates 2 by 1 matrix to do summing 
sm<-matrix(c(1,1),ncol=1)
# gets the vectors for xi from the data frame des and averages all posible pairs of points
 if (nfac>0) {x1<-des$x1; a<-combn(x1,2); x1n<-(t(a)%*%sm)/2; x1<-c(x1,x1n)} 
 if (nfac>1) {x2<-des$x2; b<-combn(x2,2); x2n<-(t(b)%*%sm)/2; x2<-c(x2,x2n)}
 if (nfac>2) {x3<-des$x3; c<-combn(x3,2); x3n<-(t(c)%*%sm)/2; x3<-c(x3,x3n)}
 if (nfac>3) {x4<-des$x4; c<-combn(x4,2); x4n<-(t(c)%*%sm)/2; x4<-c(x4,x4n)}
 if (nfac>4) {x5<-des$x5; c<-combn(x5,2); x5n<-(t(c)%*%sm)/2; x5<-c(x5,x5n)}
 if (nfac>5) {x6<-des$x6; c<-combn(x6,2); x6n<-(t(c)%*%sm)/2; x6<-c(x6,x6n)}
 if (nfac>6) {x7<-des$x7; c<-combn(x7,2); x7n<-(t(c)%*%sm)/2; x7<-c(x7,x7n)}
 if (nfac>7) {x8<-des$x8; c<-combn(x8,2); x8n<-(t(c)%*%sm)/2; x8<-c(x8,x8n)}
 if (nfac>8) {x9<-des$x9; c<-combn(x9,2); x9n<-(t(c)%*%sm)/2; x9<-c(x9,x9n)}
 if (nfac>9) {x10<-des$x10; c<-combn(x10,2); x10n<-(t(c)%*%sm)/2; x10<-c(x10,x10n)}
 if (nfac>10) {x11<-des$x11; c<-combn(x11,2); x11n<-(t(c)%*%sm)/2; x11<-c(x11,x11n)}
 if (nfac>11) {x12<-des$x12; c<-combn(x12,2); x12n<-(t(c)%*%sm)/2; x12<-c(x12,x12n)}
# makes new data frame containing original plus averages of all pairs of points
 if (nfac==2) {new<-data.frame(x1,x2)}
 if (nfac==3) {new<-data.frame(x1,x2,x3)}
 if (nfac==4) {new<-data.frame(x1,x2,x3,x4)}
 if (nfac==5) {new<-data.frame(x1,x2,x3,x4,x5)}
 if (nfac==6) {new<-data.frame(x1,x2,x3,x4,x5,x6)}
 if (nfac==7) {new<-data.frame(x1,x2,x3,x4,x5,x6,x7)}
 if (nfac==8) {new<-data.frame(x1,x2,x3,x4,x5,x6,x7,x8)}
 if (nfac==9) {new<-data.frame(x1,x2,x3,x4,x5,x6,x7,x8,x9)}
 if (nfac==10) {new<-data.frame(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)}
 if (nfac==11) {new<-data.frame(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11)}
 if (nfac==12) {new<-data.frame(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12)}
# remove duplicate rows from the new data frame
new<-new[!duplicated(new),]
return(new)
 }

