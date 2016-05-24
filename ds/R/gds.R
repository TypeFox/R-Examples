gds <-
function(data){
nam=names(as.data.frame(data)) 
x=as.data.frame(data)
d1=sapply(x,list)
d2=sapply(d1,na.omit)
h=ifelse(is.list(d2),1,2)
ww=list(f1=function(d2){return(d2)}, f2=function(d2){return(as.data.frame(d2))})
d5=ww[[h]]
d5=d5(d2)
qww=ifelse(ncol(x)==1,2,1)
qw=list(
n1=sapply(d5,length),
n2=nrow(d2))
n=qw[[qww]]
f1=function(x){r=sum((x-mean(x, na.rm=TRUE))^2, na.rm=T) ;return(r)}
s1= sapply(x,f1)
mm2=matrix(s1,ncol=1)
mm2=cbind(mm2,n)
m2=mm2[,1]/mm2[,2]
f2=function(x){r=sum((x-mean(x, na.rm=TRUE))^3, na.rm=T) ;return(r)}
s2= sapply(x,f2)
mm3=matrix(s2,ncol=1)
mm3=cbind(mm3,n)
m3=mm3[,1]/mm3[,2]
f3=function(x){r=sum((x-mean(x, na.rm=TRUE))^4, na.rm=T) ;return(r)}
s3= sapply(x,f3)
mm4=matrix(s3,ncol=1)
mm4=cbind(mm4,n)
m4=mm4[,1]/mm4[,2]
med=function(x){m=mean(x,na.rm=TRUE); return(m)}
maxi=function(x){maxi=max(x,na.rm=TRUE); r=return(maxi)}
mini=function(x){mini=min(x,na.rm=TRUE); r=return(mini)}
medi=function(x){medi=median(x,na.rm=TRUE); return(medi)}
quant=function(x){quant=quantile(x,na.rm=TRUE, prob=c(0.9987,0.9773,0.8414,0.1587,0.0228,0.0014)); return(quant)}
p=sapply(x, shapiro.test)
m=sapply(x,med)
max=sapply(x,maxi) 
min=sapply(x,mini) 
medi=sapply(x,medi) 
quantil=sapply(x,quant) 
a=max-min
vari= as.numeric(mm2[,1])/as.numeric(mm2[,2]-1)
ste=vari^0.5
mse=ste/(n^0.5)
cv=ste*100/m
as=m3/(m2^(3/2))
cur=m4/m2^2
p=as.numeric(p[2,])
med1=m+ste
med2=m+2*ste
med3=m+3*ste
med11=m-ste
med22=m-2*ste
med33=m-3*ste
r=rbind(m,max,min,medi,med3,med2,med1,med11,med22,med33,quantil,n, a,vari,ste,mse,cv,as,cur,p)
r=round(r,4)
r=as.data.frame(r)
rownames(r)=c("Mean", "Maximum","Minimum","Median","Mean + 3 standard deviation","Mean + 2 standard deviation","Mean + 1 standard deviation","Mean - 1 standard deviation","Mean - 2 standard deviation","Mean - 3 standard deviation","Quantile (99.87%)","Quantile (97.73%)","Quantile (84.14%)","Quantile (15.87%)","Quantile (2.28%)","Quantile (0.14%)","n","Range","Variance", "Standard deviation", "Standard error of the mean","Coefficient of variation (%)", "Skewness", "Kurtosis", "P-value (Shapiro-Wilk)")
colnames(r)=nam
return(r)
}
