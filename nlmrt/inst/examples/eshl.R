#,"ESHL SCHEUNEMAN  PROBLEM -- 851111"
rm(list=ls())
xx<-c(2.4,1.7,3.2,5.0,3.2,2.0,4.1,4,3,1,6.2)
yy<-c(7,15.2,10.4,21.5,16.8,24.8,21.4,5,5.2,2.9,11.2)
esdata<-data.frame(xx=xx, yy=yy)
resx<-yy~B1*(xx^B2)
st1<-c(B1=1, B2=1)
require(nlmrt)
anlsmnq1<-nlxb(resx, start=st1, trace=TRUE, data=esdata)
print(anlsmnq1)
anls1<-nls(resx, start=st1, trace=TRUE, data=esdata)
print(anls1)

