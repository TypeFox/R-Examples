forecastGARCH <-
function(fitARMA,fitGARCH,r=3,trace=FALSE,newxreg=NULL){
if(!is.null(newxreg)) 
{
if(is.vector(newxreg)) if(length(newxreg)!=1) stop("Length xreg must be 1!")
if(is.data.frame(newxreg)) if(dim(newxreg)[1]!=1) stop("Length xreg must be 1!")
}

namthang<-function(data.ts){
batdau<-start(data.ts)
tanso<-frequency(data.ts)
nam1<-batdau[1]
thang1<-batdau[2]
ketthuc<-end(data.ts)
nam2<-ketthuc[1]
thang2<-ketthuc[2]
namkq<-1:length(data.ts)
thangkq<-1:length(data.ts)
index=0;
for(nam in nam1:nam2)
for(thang in 1:tanso)
if(nam!=nam1 || thang>=thang1){
index=index+1
namkq[index]<-nam
thangkq[index]<-thang
if(nam==nam2 & thang==thang2) break;
}
if(tanso==4)
{thangkq[thangkq==1]<-"Q1";thangkq[thangkq==2]<-"Q1";
thangkq[thangkq==3]<-"Q1";thangkq[thangkq==4]<-"Q1";
print<-paste(namkq,thangkq,sep=" ")}
else
if(tanso==12)
{thangkq[thangkq==1]<-"Jan";thangkq[thangkq==2]<-"Feb";
thangkq[thangkq==3]<-"Mar";thangkq[thangkq==4]<-"Apr";
thangkq[thangkq==5]<-"May";thangkq[thangkq==6]<-"Jun";
thangkq[thangkq==7]<-"Jul";thangkq[thangkq==8]<-"Aug";
thangkq[thangkq==9]<-"Sep";thangkq[thangkq==10]<-"Oct";
thangkq[thangkq==11]<-"Nov";thangkq[thangkq==12]<-"Dec";
print<-paste(namkq,thangkq,sep=" ")}
else
if(tanso==7)
{thangkq[thangkq==1]<-"Mon";thangkq[thangkq==2]<-"Tue";
thangkq[thangkq==3]<-"Wed";thangkq[thangkq==4]<-"Thu";
thangkq[thangkq==5]<-"Fri";thangkq[thangkq==6]<-"Sat";
thangkq[thangkq==7]<-"Sun";
print<-paste(namkq,thangkq,sep=" ")}
else
if(tanso!=1)
print<-paste("(",namkq,",",thangkq,")",sep="")
else
print<-namkq
print
}

swap<-function(ts){
x<-1:length(ts)
for(i in 1:length(ts))x[i]<-ts[length(ts)-i+1]
x
}


if(class(fitARMA)!="Arima")stop("Error for fitARMA!")
if(class(fitGARCH)!="garch")stop("Error for fitGARCH!")


coefGARCH<-fitGARCH$coef
n=length(coefGARCH)
vt.garch<-0
for(i in 1:n) if(names(coefGARCH[i])=="b1") {vt.garch<-i;break}
if(vt.garch!=0){
cfARCH<-coefGARCH[1:(vt.garch-1)]
cfGARCH<-coefGARCH[vt.garch:length(coefGARCH)]
}
else
cfARCH<-coefGARCH

if(vt.garch!=0) condition<-c(vt.garch-1) else condition<-length(coefGARCH)
resARMA<-fitARMA$resid
res<-resARMA[(length(resARMA)-condition+2):length(resARMA)]
res2<-res^2
res2<-swap(res2)
temp.E<-cfARCH[-1]*res2

if(vt.garch!=0){
n.var.necessarry<-length(cfGARCH)
varGARCH<-fitGARCH$fitted.value[,1]
var<-varGARCH[(length(varGARCH)-n.var.necessarry+1):length(varGARCH)]
var2<-var^2
var2<-swap(var2)
temp.V<-cfGARCH*var2
}

if(vt.garch!=0)Var<-sum(c(cfARCH[1],temp.E,temp.V)) else Var<-sum(c(cfARCH[1],temp.E))
SSL<-as.numeric(predict(fitARMA,n.ahead=1,newxreg=newxreg)$pred)

if(vt.garch!=0)n.row<-max((length(cfARCH)-1),length(cfGARCH))else n.row<-(length(cfARCH)-1)

ts.temp<-ts(rep(1,(length(resARMA)+1)),start=start(resARMA),frequency=frequency(resARMA))
point<-namthang(ts.temp)[(length(ts.temp)-n.row):length(ts.temp)]

kq.point<-point
kq.res<-c(rep(" ",(n.row-length(res))),round(res,r)," ")
kq.sq.res<-c(rep(" ",(n.row-length(res2))),round(res2,r)," ")
kq.ssl<-c(rep(" ",n.row),round(SSL,r))
kq.var<-c(rep(" ",n.row),round(Var,r))

if(vt.garch!=0){
var2<-swap(var2)
temp<-(length(var2)+1)
for(i in n.row:(n.row-length(cfGARCH)+1))
{temp=temp-1
kq.var[i]<-round(var2[temp],r)
}}

kq<-data.frame(kq.point,kq.res,kq.sq.res,kq.ssl,kq.var)
dimnames(kq)[[2]]<-c("Point","res","res^2","SSL.forecast","VAR.forecast")
anser<-list(ARCH=coefGARCH,ARMA=fitARMA$coef,forecast=kq)
kq.little<-kq[dim(kq)[1],c(1,4,5)]
dimnames(kq.little)[[1]]<-""

if(trace==FALSE) print(kq.little) else if(trace==TRUE)print(anser) else stop("'trace' must is TRUE or FALSE!")
}
