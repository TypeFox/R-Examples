`JarqueBeraTest` <-
function(z){
x<-z-mean(z)
n<-length(x)
m2<-sum(x^2)/n
m3<-sum(x^3)/n
m4<-sum(x^4)/n
g3<-m3/(m2^(3/2))
g4<-m4/(m2^2)
LM<-n*((g3^2)/6+((g4-3)^2)/24)
pv<-1-pchisq(LM,2)
list(LM=LM,pvalue=pv)
}

