# Put in ref to Hock and Schitkowski problem 25
# 3000 PRINT "PROBLEM HS25   851001"
m<-readline("Number of data points (generally 99):")
sdevn<-readline("standard deviation of disturbances:")
if (sdevn > 0) y1<-rnorm(m, mean=0, sd=sdevn) else y1<-rep(0,m)
iseq<-seq(1,m)

ymod<-0.01*iseq
y2<-25+(-50*log(ymod))^(2/3)
formula<- ymod ~ exp((-1/b1)*((y2-b2)^b3))
lo<-c(0.1, 0, 0)
up<-c(100, 25.6, 5)
st<-c(b1=100, b2=12.5, b3=3)
# st<-c(b1=50, b2=25, b3=1.5) # solution

hdata<-data.frame(ymod=ymod, y2=y2)

ahs25nlsmnqb<-nlsmnqb(formula, st, trace=TRUE, data=hdata, lower=lo, upper=up, control=list(watch=TRUE))

print(ahs25nlsmnqb)
