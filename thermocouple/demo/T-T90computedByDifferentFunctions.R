# T-T90 as computed by different functions
# Franco Pavese and Gianfranco Molinar Min Beciet, 2013
# Modern Gas-Based Temperature and Pressure Measurements
# Springer Science
# Fig. 1.6, pp. 42

C=c(-0.003073048 
,0.002521906 
,- 0.007092070 
,- 0.007160607 
,0.011774420 
,0.015114007 
,0.000330789 
,0.060097369 
,0.042668756 )
K=c(2.3, 2.3, 2.3, 2.3, 90, 300, 300, 450, 450, 1238, 1238, 1238, 1238)
xL<-c(200, 1200)
yL<-c(-.01, .06)
a <- as.numeric(unlist(TminusT90Pavese4CubicPolynomials(200:1200)))
b <- as.numeric(unlist(TminusT90CCT2008(200:1200)))
C <- as.numeric(unlist(SplineEval(200:1200,K,C)))

plot(200:1200, C, type='l',xlim=xL, ylim=yL, col='cyan',xlab='T/K',ylab='(T-T90)/K',main='T-T90 as computed by different functions');par(new=T)
plot(200:1200, b, type='l',xlim=xL, ylim=yL, col='red',xlab='',ylab='');par(new=T)
plot(200:1200, abs(C - b), type='l', xlim=xL, ylim=yL, col='gray',xlab='',ylab='')
legend('topleft', c('Four subintervals spline','Polynomial (CCT WG4 2008)','Diff of the above'), text.col=c('cyan','red','gray'))
