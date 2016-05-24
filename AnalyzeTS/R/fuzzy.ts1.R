fuzzy.ts1 <-
function(ts,n=5,D1=0,D2=0,type=c("Chen","Singh","Heuristic","Chen-Hsu"),bin=NULL,trace=FALSE,divide=NULL,plot=FALSE){

#--ham con---
is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  
abs(x - round(x)) < tol

matdo<-function(x,g){
method <- 1
x.unique <- sort(unique(x))
min.dif <- min(diff(x.unique),na.rm=TRUE)/2
min.dif.factor <- 1
nnm <- sum(!is.na(x))
n <- table(x)
xx <- as.double(names(n))
cum <- cumsum(n)
m <- length(xx)
y <- as.integer(ifelse(is.na(x), NA, 1))
cuts <- approx(cum,xx,xout=(1:g)*nnm/g,method="constant",rule=2,f=1)$y
cuts[length(cuts)] <- max(xx)
lower <- xx[1]
upper <- 1e+45
up <- low <- double(g)
i <- 0
for (j in 1:g) {
cj <- if (method == 1 || j == 1) cuts[j]
else {
s <- if (is.na(lower)) 
FALSE
else xx >= lower
cum.used <- if (all(s)) 0
else max(cum[!s])
if (j == m) max(xx)
else if (sum(s) < 2) max(xx)
else approx(cum[s] - cum.used, xx[s], xout = (nnm - 
cum.used)/(g - j + 1), method = "constant", rule = 2, f = 1)$y
}
if (cj == upper) next
i <- i + 1
upper <- cj
y[x >= (lower - min.dif.factor * min.dif)] <- i
low[i] <- lower
lower <- if (j == g) upper
else min(xx[xx > upper],na.rm=TRUE)
if (is.na(lower)) lower <- upper
up[i] <- lower
}
low <- low[1:i]
up <- up[1:i]
kq=c(low,up[length(up)])
v=rep(0,length(kq))
for(i in 1:length(kq))
{tt=length(x[x==kq[i]])
if(tt!=0)v[i]=i}
v[v==0]<-NA
v<-na.omit(v)
for(i in 1:length(v))
kq[v[i]]=kq[v[i]]-10^(-3)
kq
}

rule1<-function(x2,x1,A3,D){
#trich thong tin
tapmo<-as.character(D[,1])
low<-D[,2]
up<-D[,3]
bw<-D[,4]
for(i in 1:length(tapmo)) #dinh vi tri A3
if(A3==tapmo[i]) {vt=i;break}
s<-(up[vt]-low[vt])/4 #tinh khoang chia
different<-abs(x2-x1)/2
distant=abs(up[vt]-low[vt])/2
if(different > distant) x=low[vt]+3*s
else
if(different == distant) x=low[vt]+2*s
else
if(different < distant) x=low[vt]+1*s
x
}

rule2<-function(x3,x2,x1,A4,D){
#trich thong tin
tapmo<-as.character(D[,1])
low<-D[,2]
up<-D[,3]
bw<-D[,4]
for(i in 1:length(tapmo)) #dinh vi tri A2
if(A4==tapmo[i]) {vt=i;break}
s<-(up[vt]-low[vt])/4 #tinh khoang chia
different1a<- abs((x3-x2)-(x2-x1))*2+x3
different1b<- x3-abs((x3-x2)-(x2-x1))*2
different2a<- abs((x3-x2)-(x2-x1))/2+x3
different2b<- x3-abs((x3-x2)-(x2-x1))/2
if(low[vt]<=different1a &  different1a<=up[vt]) x=low[vt]+3*s
else
if(low[vt]<=different1b &  different1b<=up[vt]) x=low[vt]+3*s
else
if(low[vt]<=different2a &  different2a<=up[vt]) x=low[vt]+1*s
else
if(low[vt]<=different2b &  different2b<=up[vt]) x=low[vt]+1*s
else
x=low[vt]+2*s
x
}

rule3<-function(x3,x2,x1,A4,D){
#trich thong tin
tapmo<-as.character(D[,1])
low<-D[,2]
up<-D[,3]
bw<-D[,4]
for(i in 1:length(tapmo)) #dinh vi tri A2
if(A4==tapmo[i]) {vt=i;break}
s<-(up[vt]-low[vt])/4 #tinh khoang chia
different1a<- abs((x3-x2)-(x2-x1))*2+x3
different1b<- x3-abs((x3-x2)-(x2-x1))*2
different2a<- abs((x3-x2)-(x2-x1))/2+x3
different2b<- x3-abs((x3-x2)-(x2-x1))/2
if(low[vt]<=different2a &  different2a<=up[vt]) x=low[vt]+1*s
else
if(low[vt]<=different2b &  different2b<=up[vt]) x=low[vt]+1*s
else
if(low[vt]<=different1a &  different1a<=up[vt]) x=low[vt]+3*s
else
if(low[vt]<=different1b &  different1b<=up[vt]) x=low[vt]+3*s
else
x=low[vt]+2*s
x
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

divide.distance<-function(tsi,f,min.tsi,max.tsi){
min.x=min(tsi,min.tsi)
max.x=max(tsi,max.tsi)
h=(max.x-min.x)/f
k<-1:(f+1);
for(i in 1:(f+1)){
if(i==1)k[i]=min.x
else
k[i]=min.x+(i-1)*h
}
k
}

quanhe<-function(ts,type){
qh<-1:length(ts)
if(type=="Chen") for(i in 1:length(ts))if(i>1)qh[i]<-c("<--")
if(type=="Singh") for(i in 1:length(ts))if(i>1)qh[i]<-c("<--")
if(type=="Heuristic") for(i in 1:length(ts))if(i>1)qh[i]<-c("<--")
if(type=="Chen-Hsu") for(i in 1:length(ts))if(i>1)qh[i]<-c("<--")
qh[qh!="<--"]<-"-x-"
qh}
#--ham con---


#---kiem tra input---
#Kiem tra so lieu
if (!is.numeric(ts)) stop("Error in 'ts'!")

if(!is.ts(ts)) stop("Error in 'ts'!") else if(!is.null(dim(ts)))stop("Error in 'ts'!")

kt<-0
for(i in 1:length(ts))if(is.na(ts[i]))kt=kt+1
if(kt>0)stop("Trong chuoi co gia tri NA!")

if (is.na(D1)||!is.numeric(D1)) stop("Error in 'D1'!")
if (is.na(D2)||!is.numeric(D2)) stop("Error in 'D2'!")
if (is.na(n)||!is.numeric(n) || n < 1 || !is.wholenumber(n)) stop("Error in 'n'!")

type<-match.arg(type)

if(type!="Chen" & type!="Singh" & type!="Heuristic" & type!="Chen-Hsu") stop("Error in 'type'!")

if(type!="Chen-Hsu")if(!is.null(bin))stop("'bin' just for Chen-Hsu model!")
if(!is.null(divide))if(type!="Chen-Hsu") stop("'divide' just for Chen-Hsu model!")
if(type=="Chen-Hsu")if(is.null(divide))divide<-c("distance") else if(divide!="distance" & divide!="density")stop("'divide' must be 'distance' or 'density'!")
#---kiem tra input---


#---tinh toan ban dau---
d1=D1
d2=D2
thoidiem<-namthang(ts)
quanhemo<-quanhe(ts,type)

#Buoc 1: Xay dung va chia tap nen
min.x=min(ts)-d1
max.x=max(ts)+d2
h=(max.x-min.x)/n
k<-1:(n+1);U<-1:n
for(i in 1:(n+1)){
if(i==1)k[i]=min.x
else
{k[i]=min.x+(i-1)*h
 U[i-1]=paste("A",i-1,sep="")}
}

D<-data.frame(U,low=k[1:n],up=k[2:(n+1)])
D$Bw<-(1/2)*(D$low+D$up)

#Buoc 2: Xac dinh cac tap mo
loai<-1:length(ts)
for(i in 1:length(ts)){
for(j in 1:n){
if (D$low[j]<=ts[i] & ts[i]<=D$up[j])
{loai[i]=paste("A",j,sep=""); break;}
}}

#Buoc 3: Xac dinh moi quan he mo
loai.old<-1:length(loai)
for(i in 1:length(loai)){
if(i==1)loai.old[i]=NA
else
loai.old[i]=loai[i-1]
}
quanhemokq<-paste(loai,quanhemo,loai.old,sep="")
D1<-data.frame(ts,loai,loai.old)
b<-table(D1$loai)
ni<-1:n
for(i in 1:n)
for(j in 1:n){
if(paste("A",i,sep="")==names(b)[j] & j<=length(b)) 
{ni[i]=b[j];break;}else{
if(j>length(b)) {ni[i]=0; break;}
}
}
D$ni<-ni

#Gom nhom
a<-1:(n*n)
a<-matrix(a,nrow=n)

for(cot in 1:n){
for(i in 2:length(D1$loai.old)){
for(j in 1:n){
if(D1$loai.old[i]==paste("A",cot,sep="") & a[j,cot]!=paste("A",j,sep="")
    & D1$loai[i]==paste("A",j,sep="")) {a[j,cot]=paste("A",j,sep=""); break;}
}}}
for(cot in 1:n){
for(j in 1:n) if(a[j,cot]!=paste("A",j,sep="")) a[j,cot]<-NA
}
#---tinh toan ban dau---


#---giai mo va du bao---

#type Chen
if(type=="Chen"){
#1.Giai mo
mo<-1:n
for(cot in 1:n){
if(cot>1)a[a==1]<-NA
s<-length(na.omit(a[,cot]))
a[is.na(a)]<-1
if(s==1){
for(j in 1:n) if(a[j,cot]==paste("A",j,sep="")){t=D$Bw[j];break;}
mo[cot]=t
}
else
if(s>1){
t=0
for(j in 1:n) if(a[j,cot]== paste("A",j,sep=""))t=t+D$Bw[j]
mo[cot]=t/s
}
else
if(s==0){
mo[cot]=D$Bw[cot]
}
}

#2.Du bao
db<-1:length(loai.old)
for(i in 1:length(loai.old)){
if(i==1)db[i]=NA
else
for(j in 1:n)if(loai.old[i]==paste("A",j,sep="")) {db[i]=mo[j];break;}
}
db<-ts(db,start=start(ts),frequency=frequency(ts))
DB1<-data.frame(thoidiem,D1$ts,quanhemokq,db)
}#finish Chen

#type == Singh
if(type=="Singh"){
F<-1:length(ts);Ds<-1:length(ts)
X<-1:length(ts);XX<-1:length(ts)
Y<-1:length(ts);YY<-1:length(ts)
P<-1:length(ts);PP<-1:length(ts)
Q<-1:length(ts);QQ<-1:length(ts)
G<-1:length(ts);GG<-1:length(ts)
H<-1:length(ts);HH<-1:length(ts)
E<-ts;

#Xoa bo level thua: begin
tt=0;v<-1:length(D$ni)
for(t0 in 1:length(D$ni))
if(D$ni[t0]==0){tt=tt+1;v[tt]=t0;}
if(tt==0)(v=NULL) else v=v[1:tt]

if(!is.null(v)){
Dt<-D[-v,]
t1=1;
while(t1<=length(levels(Dt$U))){
test=0;
for(t2 in 1:length(Dt$U))
if(levels(Dt$U)[t1]==Dt$U[t2]){test=1;break;}

if(test==0)levels(Dt$U)[t1]=NA else t1=t1+1;
}}

if(is.null(v))Dt<-D
#Xoa bo level thua: the end

for(i in 0:(length(ts)-1)){
if(i<3)F[i+1]<-NA
else{
R=0;S=0;
Ds[i]=abs(abs(E[i]-E[i-1])-abs(E[i-1]-E[i-2]))
X[i]=E[i]+(Ds[i]/2)
XX[i]=E[i]-(Ds[i]/2)
Y[i]=E[i]+Ds[i]
YY[i]=E[i]-Ds[i]
P[i]=E[i]+(Ds[i]/4)
PP[i]=E[i]-(Ds[i]/4)
Q[i]=E[i]+2*Ds[i]
QQ[i]=E[i]-2*Ds[i]
G[i]=E[i]+(Ds[i]/6)
GG[i]=E[i]-(Ds[i]/6)
H[i]=E[i]+3*Ds[i]
HH[i]=E[i]-3*Ds[i]

Dt.U<-as.character(Dt$U)
D1.loai<-as.character(D1$loai)

for(j in 1:length(Dt.U)) if(D1.loai[i+1]==Dt.U[j]) {l=j;break;}

if(Dt$low[l]<=X[i] & X[i]<=Dt$up[l]) {R=R+X[i];S=S+1;}
if(Dt$low[l]<=XX[i] & XX[i]<=Dt$up[l]) {R=R+XX[i];S=S+1;}

if(Dt$low[l]<=Y[i] & Y[i]<=Dt$up[l]) {R=R+Y[i];S=S+1;}
if(Dt$low[l]<=YY[i] & YY[i]<=Dt$up[l]) {R=R+YY[i];S=S+1;}

if(Dt$low[l]<=P[i] & P[i]<=Dt$up[l]) {R=R+P[i];S=S+1;}
if(Dt$low[l]<=PP[i] & PP[i]<=Dt$up[l]) {R=R+PP[i];S=S+1;}

if(Dt$low[l]<=Q[i] & Q[i]<=Dt$up[l]) {R=R+Q[i];S=S+1;}
if(Dt$low[l]<=QQ[i] & QQ[i]<=Dt$up[l]) {R=R+QQ[i];S=S+1;}

if(Dt$low[l]<=G[i] & G[i]<=Dt$up[l]) {R=R+G[i];S=S+1;}
if(Dt$low[l]<=GG[i] & GG[i]<=Dt$up[l]) {R=R+GG[i];S=S+1;}

if(Dt$low[l]<=H[i] & H[i]<=Dt$up[l]) {R=R+H[i];S=S+1;}
if(Dt$low[l]<=HH[i] & HH[i]<=Dt$up[l]) {R=R+HH[i];S=S+1;}

F[i+1]=(R + Dt$Bw[l])/(S+1)
}
}
F<-ts(F,start=start(ts),frequency=frequency(ts))
DB2<-data.frame(thoidiem,sl.goc=D1$ts,quanhemo=quanhemokq,.....sl.mo=F)
}#finish Singh

#type Heuristic
if(type=="Heuristic"){
#Xac dinh nhom quan he mo Heuristic 
tap="begin";t.thai="calm";Bw.h=0

for(cot in 1:n){
locate<-rep(0,n)
for(i in 1:n) if(!is.na(a[i,cot]))locate[i]=i
locate[locate==0]<-NA
locate<-na.omit(locate)

if(length(locate)>0){
#xu huong giam: lon->nho
l=0
locate[length(locate)+1]<-(n+1) 
for(i in 1:(length(locate)-1)) 
  if(cot >= locate[i] & cot < locate[i+1]) 
     {l=i;break;}
locate=locate[1:(length(locate)-1)]
if(l!=0){
t=0  
for(i in 1:l) t=t+D$Bw[locate[i]]
tb<-t/length(locate[1:l])
tap[length(tap)+1]=paste("A",cot,sep="")
t.thai[length(t.thai)+1]="down"
Bw.h[length(Bw.h)+1]=tb
}

#xu huong tang: nho->lon
l=0
for(i in 1:length(locate)) if(cot <= locate[i]) {l=i;break;}
if(l!=0){
t=0  
for(i in l:length(locate)) t=t+D$Bw[locate[i]]
tb<-t/length(locate[l:length(locate)])
tap[length(tap)+1]=paste("A",cot,sep="")
t.thai[length(t.thai)+1]="raise"
Bw.h[length(Bw.h)+1]=tb
}}}
tap<-factor(tap[2:length(tap)])
t.thai<-t.thai[2:length(t.thai)]
Bw.h<-Bw.h[2:length(Bw.h)]
D.h<-data.frame(tap,t.thai,Bw.h)

#Du bao
Heuristic<-1:length(ts)
Heuristic[1]<-NA

for(j in 2:length(ts)){ 
xj<-ts[j]-ts[j-1]

D1.loai.old<-as.character(D1$loai.old)
D.h.tap<-as.character(D.h$tap)

if(xj>=0){
for(v in 1:dim(D.h)[1])
if(D1.loai.old[j]==D.h.tap[v] & D.h$t.thai[v]=="raise") bt=v
Heuristic[j]=D.h$Bw.h[bt]
}
else
if(xj<0){
for(v in 1:dim(D.h)[1])
if(D1.loai.old[j]==D.h.tap[v] & D.h$t.thai[v]=="down") bt=v
Heuristic[j]=D.h$Bw.h[bt]
}
}
Heuristic<-ts(Heuristic,start=start(ts),frequency=frequency(ts))
DB3<-data.frame(thoidiem,sl.goc=D1$ts,quanhemo=quanhemokq,.....sl.mo=Heuristic)
}#finish Heuristic

#type Chen va Hsu
if(type=="Chen-Hsu"){
if(is.null(bin)){ #bin==null
#Phan lai tap mo
b<-D$ni
level=max(round(length(ts)/n-0.5),min(b),2)#tuy theo y nguoi dung
f<-1:n
for(i in 1:n) {if(b[i]>(level+2))f[i]=round(b[i]/level+0.5) else f[i]=1;
               if(f[i]==0) f[i]=1}
n2=sum(f)
U2<-1:n2; k2<-1:(n2+1); t=0; k2[1]<-k[1]
for( i in 1:n){
tsi=ts[ts>=D$low[i]]
tsi=tsi[tsi<=D$up[i]]
if(divide=="density")new=matdo(tsi,f[i])
else
if(divide=="distance")new=divide.distance(tsi,f[i],D$low[i],D$up[i])
else
stop("divide phai co gia tri la 'distance' hoac 'frequency'!")
new[1]=min(new[1],D$low[i])
new[length(new)]=max(new[length(new)],D$up[i])
for(j in 2:(f[i]+1)){
t=t+1
k2[t+1]=new[j]
U2[t]=paste("A",t,sep="")
}
}
} else{   #bin !=null
if(min(bin)<=min.x || max(bin)>=max.x || length(bin)<1)
stop("gia tri tam bin khong hop le!")
k2<-1:(length(bin)+2);U2<-1:(length(bin)+1)
k2[1]<-min.x;k2[length(k2)]<-max.x
for(t in 2:(length(k2)-1)) k2[t]=bin[t-1]
for(t in 1:(length(bin)+1)) U2[t]=paste("A",t,sep="")
n2<-length(U2)
}

D.ch<-data.frame(U2,low2=k2[1:n2],up2=k2[2:(n2+1)])
D.ch$Bw2<-(1/2)*(D.ch$low+D.ch$up)

#Xac dinh moi quan he mo moi
loai2<-1:length(ts)
for(i in 1:length(ts)){
for(j in 1:n2){
if (D.ch$low2[j]<=ts[i] & ts[i]<=D.ch$up2[j])
{loai2[i]=paste("A",j,sep=""); break;}
}}
D1.ch<-data.frame(ts,loai2,loai2.old=c(NA,loai2[1:(length(loai2)-1)]))
b<-table(D1.ch$loai2)
quanhemokq<-paste(D1.ch$loai2,quanhemo,D1.ch$loai2.old,sep="")

ni<-1:n2
for(i in 1:n2)
for(j in 1:n2){
if(paste("A",i,sep="")==names(b)[j] & j<=length(b)) 
{ni[i]=b[j]; break;}else{
if(j>length(b)){ni[i]=0;break;}
}
}
D.ch$ni<-ni

#Xoa bo level thua: begin
tt=0;v<-1:length(D.ch$ni)
for(t0 in 1:length(D.ch$ni))
if(D.ch$ni[t0]==0){tt=tt+1;v[tt]=t0;}
if(tt==0)(v=NULL) else v=v[1:tt]
if(!is.null(v)){
Dt<-D.ch[-v,]
t1=1;
while(t1<=length(levels(Dt$U2))){
test=0;
for(t2 in 1:length(Dt$U2))
if(levels(Dt$U2)[t1]==Dt$U2[t2]){test=1;break;}
if(test==0)levels(Dt$U2)[t1]=NA else t1=t1+1;
}}
else
Dt<-D.ch
#Xoa bo level thua: the end

#Giai mo va du bao
X=ts
A<-as.character(D1.ch$loai2)
P<-1:length(X)

tapmochonam2<-as.character(Dt[,1])
for(i in 1:length(tapmochonam2)) #dinh vi tri tap mo cua nam 2
if(A[2]==tapmochonam2[i]) {vt=i;break} 
P[1]<-NA;P[2]=Dt[vt,4]

for(t in 2:(length(ts)-1)){
if(t==2) P[t+1]=rule1(X[t],X[t-1],A[t+1],Dt)
else
{
different=(X[t]-X[t-1])-(X[t-1]-X[t-2])
if(A[t+1]>A[t] & different>=0) P[t+1]=rule2(X[t],X[t-1],X[t-2],A[t+1],Dt)#th1
else
if(A[t+1]>A[t] & different<=0) P[t+1]=rule3(X[t],X[t-1],X[t-2],A[t+1],Dt)#th2
else
if(A[t+1]<A[t] & different>=0) P[t+1]=rule2(X[t],X[t-1],X[t-2],A[t+1],Dt)#th3
else
if(A[t+1]<A[t] & different<=0) P[t+1]=rule3(X[t],X[t-1],X[t-2],A[t+1],Dt)#th4
else
if(A[t+1]==A[t] & different>=0) P[t+1]=rule2(X[t],X[t-1],X[t-2],A[t+1],Dt)#th5
else
if(A[t+1]==A[t] & different<=0) P[t+1]=rule3(X[t],X[t-1],X[t-2],A[t+1],Dt)#th6
}
}

P<-ts(P,start=start(ts),frequency=frequency(ts))
DB4<-data.frame(thoidiem,sl.goc=D1.ch$ts,quanhemo=quanhemokq,".....sl.mo"=P)

diem25<-1:length(D.ch[,3])
diem75<-1:length(D.ch[,3])
h<-(D.ch[,3]-D.ch[,2])/4
canduoi<-D.ch[,2]
for(i in 1:length(D.ch[,3])){
diem25[i]<-canduoi[i]+1*h[i]
diem75[i]<-canduoi[i]+3*h[i]}

D.ch<-data.frame(D.ch[,1],D.ch[,2],diem25,D.ch[,4],diem75,D.ch[,3],D.ch[,5])

namescot1<-"set"
namescot2<-"dow"
namescot3<-"0.25"
namescot4<-"0.50"
namescot5<-"0.75"
namescot6<-"up"
namescot7<-"num"

colnames(D.ch)<-c(namescot1,namescot2,namescot3,namescot4,namescot5,namescot6,namescot7)
}#finish chen-hsu
#---giai mo va du bao---


#---tinh toan do chinh xac---
if(type=="Chen")  accuracy<-av.res(Y=data.frame(DB1[,2]),F=data.frame(Chen=DB1[,4])) 
if(type=="Singh")  accuracy<-av.res(Y=data.frame(DB2[,2]),F=data.frame(Singh=DB2[,4]))
if(type=="Heuristic")  accuracy<-av.res(Y=data.frame(DB3[,2]),F=data.frame(Heuristic=DB3[,4]))
if(type=="Chen-Hsu") accuracy<-av.res(Y=data.frame(DB4[,2]),F=data.frame("Chen-Hsu"=DB4[,4]))
#---tinh toan do chinh xac---


#---xuat ket qua---
if(trace==TRUE){
namescot1<-"set"
namescot2<-"dow"
namescot3<-"up"
namescot4<-"mid"
namescot5<-"num"
colnames(D)<-c(namescot1,namescot2,namescot3,namescot4,namescot5)

namescot1<-"point"
namescot2<-"ts"
namescot3<-"relative"
namescot4<-"forecast"

if(type=="Chen")colnames(DB1)<-c(namescot1,namescot2,namescot3,namescot4)
if(type=="Singh")colnames(DB2)<-c(namescot1,namescot2,namescot3,namescot4)
if(type=="Heuristic")colnames(DB3)<-c(namescot1,namescot2,namescot3,namescot4)
if(type=="Chen-Hsu")colnames(DB4)<-c(namescot1,namescot2,namescot3,namescot4)

if(type=="Chen") MO<-list(type="Chen",table1=D,table2=DB1,accuracy=accuracy)
if(type=="Singh") MO<-list(type="Singh",table1=D,table2=DB2,accuracy=accuracy)
if(type=="Heuristic") MO<-list(type="Heuristic",table1=D,table2=DB3,accuracy=accuracy)
if(type=="Chen-Hsu") {
if(is.null(bin)) MO<-list(type="Chen-Hsu",table1=D,table2=D.ch,table3=DB4,accuracy=accuracy)
else 
if(!is.null(bin)) MO<-list(type="Chen-Hsu",table1=D.ch,table2=DB4,accuracy=accuracy)
}}
else
if(trace==FALSE){
if(type=="Chen") MO<-DB1[,4]
if(type=="Singh") MO<-DB2[,4]
if(type=="Heuristic") MO<-DB3[,4]
if(type=="Chen-Hsu") MO<-DB4[,4]}
else
MO<-c("trace must be 'TRUE' or 'FALSE'")
#---xuat ket qua---


#---ve bieu do du bao---
if(plot==TRUE){
if(type=="Chen") {
goc<-DB1[,2]
dubao<-DB1[,4]

if(length(goc)<50){
plot(goc,col="blue",main=paste("Chen",n,"fuzzy set"),type="o",
     pch=16,ylim=c(min(c(goc,dubao),na.rm=1),max(c(goc,dubao),na.rm=1)),
     xlab="index",ylab="data");
lines(dubao,col="red",type="o",pch=18)
legend("bottomright","(x,y)",c("ts","forecast"),col=c("blue","red"),lty=c(1,1),pch=c(16,18))
}

if(length(goc)>49){
plot(goc,col="blue",main=paste("Chen",n,"fuzzy set"),type="l",
     ylim=c(min(c(goc,dubao),na.rm=1),max(c(goc,dubao),na.rm=1)),
     xlab="index",ylab="data");
lines(dubao,col="red",type="l")
legend("bottomright","(x,y)",c("ts","forecast"),col=c("blue","red"),lty=c(1,1))
}

k.ve<-par()$yaxp
h=(k.ve[2]-k.ve[1])/k.ve[3]
for(i in 0:k.ve[3]) abline(h=k.ve[1]+h*i,lty=3,col="gray")
}


if(type=="Singh"){
goc<-DB2[,2]
dubao<-DB2[,4]

if(length(goc)<50){
plot(goc,col="blue",main=paste("Singh",n,"fuzzy set"),type="o",
     pch=16,ylim=c(min(c(goc,dubao),na.rm=1),max(c(goc,dubao),na.rm=1)),
     xlab="index",ylab="data");
lines(dubao,col="red",type="o",pch=18)
legend("bottomright","(x,y)",c("ts","forecast"),col=c("blue","red"),lty=c(1,1),pch=c(16,18))
}

if(length(goc)>49){
plot(goc,col="blue",main=paste("Singh",n,"fuzzy set"),type="l",
     ylim=c(min(c(goc,dubao),na.rm=1),max(c(goc,dubao),na.rm=1)),
     xlab="index",ylab="data");
lines(dubao,col="red",type="l")
legend("bottomright","(x,y)",c("ts","forecast"),col=c("blue","red"),lty=c(1,1))
}

k.ve<-par()$yaxp
h=(k.ve[2]-k.ve[1])/k.ve[3]
for(i in 0:k.ve[3]) abline(h=k.ve[1]+h*i,lty=3,col="gray")
}


if(type=="Heuristic") {
goc<-DB3[,2]
dubao<-DB3[,4]

if(length(goc)<50){
plot(goc,col="blue",main=paste("Heuristic",n,"fuzzy set"),type="o",
     pch=16,ylim=c(min(c(goc,dubao),na.rm=1),max(c(goc,dubao),na.rm=1)),
     xlab="index",ylab="data");
lines(dubao,col="red",type="o",pch=18)
legend("bottomright","(x,y)",c("ts","forecast"),col=c("blue","red"),lty=c(1,1),pch=c(16,18))
}

if(length(goc)>49){
plot(goc,col="blue",main=paste("Heuristic",n,"fuzzy set"),type="l",
     ylim=c(min(c(goc,dubao),na.rm=1),max(c(goc,dubao),na.rm=1)),
     xlab="index",ylab="data");
lines(dubao,col="red",type="l")
legend("bottomright","(x,y)",c("ts","forecast"),col=c("blue","red"),lty=c(1,1))
}

k.ve<-par()$yaxp
h=(k.ve[2]-k.ve[1])/k.ve[3]
for(i in 0:k.ve[3]) abline(h=k.ve[1]+h*i,lty=3,col="gray")
}


if(type=="Chen-Hsu") {
goc<-DB4[,2]
dubao<-DB4[,4]

if(length(goc)<50){
plot(goc,col="blue",main=paste("Chen-Hsu"),type="o",
     pch=16,ylim=c(min(c(goc,dubao),na.rm=1),max(c(goc,dubao),na.rm=1)),
     xlab="index",ylab="data");
lines(dubao,col="red",type="o",pch=18)
legend("bottomright","(x,y)",c("ts","forecast"),col=c("blue","red"),lty=c(1,1),pch=c(16,18))
}

if(length(goc)>49){
plot(goc,col="blue",main=paste("Chen-Hsu"),type="l",
     ylim=c(min(c(goc,dubao),na.rm=1),max(c(goc,dubao),na.rm=1)),
     xlab="index",ylab="data");
lines(dubao,col="red",type="l")
legend("bottomright","(x,y)",c("ts","forecast"),col=c("blue","red"),lty=c(1,1))
}

k.ve<-par()$yaxp
h=(k.ve[2]-k.ve[1])/k.ve[3]
for(i in -2:(k.ve[3]+2)) abline(h=k.ve[1]+h*i,lty=3,col="gray")
}
}
#---ve bieu do du bao---


MO

}
