# fungi.R
# "FUNGI.RES 3-PARAMETER LOGISTIC FUNCTION SETUP -  880412"
# 6 datasets available -- choose 1,2,3,4,5 or 6 
rm(list=ls())
a1x<-c(0, 0.1, 1, 10, 50, 100, 1000, 0, 0.1, 1, 10, 50, 100, 
       1000, 0, 0.1, 1, 10, 50, 100, 1000)
a1y<-c(9.0, 8.0, 9.0, 9.5, 11.0, 10.0, 4.0, 9.0, 10.0, 8.0,
        9.5, 10.0, 9.0, 4.0, 8.0, 9.0, 8.0, 11.0, 10.0, 9.5, 4.0)
a1data<-data.frame(xx=a1x, yy=a1y)

a2x<-c(0, 0.1, 1, 10, 50, 100, 1000, 0, 0.1, 1, 10, 50, 100,  1000,
        0, 0.1, 1, 10, 50, 100, 1000)
a2y<-c( 7.5, 10.5, 10.5, 11.5, 10.5, 9.0, 3.5, 7.0, 10.5,  10.5,
        10.0, 10.0, 9.0, 3.0, 10.5, 11.0, 11.0, 11.5, 10.5,  9.0, 4.0)
a2data<-data.frame(xx=a2x, yy=a2y)
b1x<-c(0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.2, 1.0, 10,  100,
        1000, 0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.2, 1.0,  10, 
        100, 1000, 0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.2,  1.0, 10,
        100, 1000)
b1y<-c(14.0, 13.5, 14.0, 12.5, 6.5, 6.5, 2.5, 1, 0, 0, 0,
        0, 14.0, 14.0, 14.0, 12.5, 6.0, 5.5, 0, 1.3, 0, 0,
        0 , 0, 13.5, 13.5, 14.5, 11, 6.5, 5, 2, 1, 0, 0, 0, 0)
b1data<-data.frame(xx=b1x, yy=b1y)
b2x<-c(0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.2, 1.0, 10,  100,
        1000, 0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.2, 1.0,  10,
        100, 1000, 0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.2,  1.0,
        10, 100, 1000)
b2y<-c(11.0, 11.5, 10.0, 9.5, 5.5, 5.0, 2.0, 0.5, 0, 0,
        0, 0, 11.0, 10.5, 11.0, 10.0, 5.5, 0, 0, 0.8, 0,
        0, 0, 0, 12.5, 11.5, 10.0, 9.5, 6.5, 5.0, 1.0, 0.3,
        0, 0, 0, 0)
b2data<-data.frame(xx=b2x, yy=b2y)
c1x<-c(0, 0.1, 1, 10, 50, 100, 1000, 0, 0.1, 1, 10, 50,
        100, 1000, 0, 0.1, 1, 10, 50, 100, 1000)
c1y<-c(10.5, 10.5, 10.5, 10.5, 9.0, 8.0, 4.0, 11.0,
        10.5, 10.0, 11.0, 9.5, 8.5, 4.0, 11.5, 10.5, 9.5, 10.0,
        9.0, 9.0, 4.5)
c1data<-data.frame(xx=c1x, yy=c1y)
c2x<-c(0, .01, .1, 1, 10, 100, 10000, .01, .1, 1,
        10, 100, 10000, .01, .1, 1, 10, 100, 1000)
c2y<-c(11.5, 11.5, 11.0, 12.0, 13.0, 13.0, 8.0, 11.0,
        11.5, 12.0, 11.0, 12.0, 11.5, 8.0, 11.0, 12.0,
        11.5, 12.0, 12.0) #, 11.5, 8.5)
c2data<-data.frame(xx=c2x, yy=c2y)

plot(a1data$xx, a1data$yy)
title("a1")
tmp<-readline("next")

X11()
plot(a2data$xx, a2data$yy)
title("a2")
tmp<-readline("next")

X11()
plot(b1data$xx, b1data$yy)
title("b1")
tmp<-readline("next")

X11()
plot(b2data$xx, b2data$yy)
title("b2")
tmp<-readline("next")

X11()
plot(c1data$xx, c1data$yy)
title("c1")
tmp<-readline("next")

X11()
plot(c2data$xx, c2data$yy)
title("c2")
tmp<-readline("next")


resx<-"yy~b1/(1+exp(-b2*(xx-b3)))"

st1<-c(b1=1, b2=1, b3=1)
cat("Data set a1\n")
a1nlsmnqb<-nlsmnqb(resx, start=st1, trace=TRUE, data=a1data, lower=0)
print(a1nlsmnqb)
tmp<-readline("next")

cat("Data set a2\n")
a2nlsmnqb<-nlsmnqb(resx, start=st1, trace=TRUE, data=a2data, lower=0)
print(a2nlsmnqb)
tmp<-readline("next")
cat("Data set b1\n")
b1nlsmnqb<-nlsmnqb(resx, start=st1, trace=TRUE, data=b1data, lower=0)
print(b1nlsmnqb)
tmp<-readline("next")
cat("Data set b2\n")
b2nlsmnqb<-nlsmnqb(resx, start=st1, trace=TRUE, data=b2data, lower=0)
print(b2nlsmnqb)
tmp<-readline("next")
cat("Data set c1\n")
c1nlsmnqb<-nlsmnqb(resx, start=st1, trace=TRUE, data=c1data, lower=0)
print(c1nlsmnqb)
tmp<-readline("next")
cat("Data set c2\n")
c2nlsmnqb<-nlsmnqb(resx, start=st1, trace=TRUE, data=c2data, lower=0)
print(c2nlsmnqb)
tmp<-readline("next")

