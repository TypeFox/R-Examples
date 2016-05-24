base.ts <-
function(ts,weight=NULL,type=c("mean.period","mean.point1","mean.point2",
"abs.period","abs.root","abs.mean","develop.period","develop.root","develop.mean",
"growth.period","growth.root","growth.mean"   )){
if(!is.vector(ts) & !is.ts(ts))stop("'ts' phai la mot vector so hoc chuoi thoi gian mot chieu")
if(!is.numeric(ts))stop("'ts' phai la chuoi so!")
if(!is.null(weight))
if(type!="mean.point2")stop("'weight' chi danh cho 'type=mean.point2'!")
if(sum(1*is.na(ts))!=0)stop("Trong chuoi co gia tri NA!")

#---muc do trung binh cua day so thoi ky---
if(type=="mean.period"){
kq<-sum(ts)/length(ts)
}
#---muc do trung binh cua day so thoi ky---


#---muc do trung binh cua day so thoi diem bang nhau---
if(type=="mean.point1"){
n<-length(ts)
s1<-sum(ts[-c(1,n)])
s2<-s1+1/2*(ts[1]+ts[n])
kq<-s2/(n-1)
}
#---muc do trung binh cua day so thoi diem bang nhau---


#---muc do trung binh cua day so thoi diem khong bang nhau---
if(type=="mean.point2"){
if(!is.vector(weight) & !is.ts(weight))stop("'weight' phai la mot vector so hoc chuoi thoi gian mot chieu")
if(length(weight)!=length(ts))stop("Chieu dai 'ts' va 'weight' phai bang nhau!")
kq<-sum(ts*weight)/sum(weight)
}
#---muc do trung binh cua day so thoi diem khong bang nhau---


#---so tuyet doi tung ky---
if(type=="abs.period"){
kq<-abs(diff(ts))
kq<-c(NA,kq)}
#---so tuyet doi tung ky---


#---so tuyet doi dinh goc---
if(type=="abs.root"){
kq<-rep(NA,length(ts))
for(i in 2:length(ts))
kq[i]<-ts[i]-ts[1]}
#---so tuyet doi dinh goc---


#---so tuyet doi trung binh---
if(type=="abs.mean"){
kq<-(ts[length(ts)]-ts[1])/(length(ts)-1)
}
#---so tuyet doi trung binh---


#---toc do phat trien tung ky---
if(type=="develop.period"){
kq<-rep(NA,length(ts))
for(i in 2:length(ts))
kq[i]<-ts[i]/ts[i-1]}
#---toc do phat trien tung ky---


#---toc do phat trien dinh goc---
if(type=="develop.root"){
kq<-rep(NA,length(ts))
for(i in 2:length(ts))
kq[i]<-ts[i]/ts[1]}
#---toc do phat trien dinh goc---


#---toc do phat trien trung binh---
if(type=="develop.mean"){
kq<-(ts[length(ts)]/ts[1])^(1/(length(ts)-1))
}
#---toc do phat trien trung binh---


#---toc do tang truong tung ky---
if(type=="growth.period"){
kq<-rep(NA,length(ts))
for(i in 2:length(ts))
kq[i]<-(ts[i]/ts[i-1])-1}
#---toc do tang truong tung ky---


#---toc do tang truong dinh goc---
if(type=="growth.root"){
kq<-rep(NA,length(ts))
for(i in 2:length(ts))
kq[i]<-(ts[i]/ts[1])-1}
#---toc do tang truong dinh goc---


#---toc do tang truong trung binh---
if(type=="growth.mean"){
kq<-(ts[length(ts)]/ts[1])^(1/(length(ts)-1))-1
}
#---toc do tang truong trung binh---




kq
}
