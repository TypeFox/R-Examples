`correl` <-
function(x,y,method = "pearson", alternative = "two.sided"){
n<-length(x)
if(method=="kendall"){
corr<-kendall(x,y)
stat<-corr$stat
rho<-corr$tau
if(alternative == "two.sided" ) pvalue<-corr$pvalue
if(alternative == "less" ) pvalue<-1-corr$pvalue/2
if(alternative == "greater") pvalue<-corr$pvalue/2
}
if(method=="spearman" ){
a<-rank(x)
b<-rank(y)
x<-a
y<-b
}
if ((method =="pearson") | (method=="spearman")) {
sumx<-sum(x^2)-sum(x)^2/n
sumy<-sum(y^2)-sum(y)^2/n
sumxy<-sum(x*y)-sum(x)*sum(y)/n
rho<-sumxy/sqrt(sumx*sumy)
gl<-n-2
stat<-rho*sqrt(gl)/(sqrt(1-rho^2))
if(alternative == "two.sided" ) pvalue<-2*(1-pt(abs(stat),gl))
if(alternative == "less" ) pvalue<-pt(abs(stat),gl)
if(alternative == "greater") pvalue<-1-pt(abs(stat),gl)
}
list(stat=stat,rho=rho,pvalue=pvalue)
}

