CMA <-
function(x,n=5){
L=n
#Kiem tra so lieu
if (!is.numeric(x)) stop("Error in series!")
is.wholenumber <-function(num, tol = .Machine$double.eps^0.5)  
abs(num - round(num)) < tol

if (is.na(L)||!is.numeric(L)||is.null(L)||!is.wholenumber(n)||L < 2) 
  stop("Error in n!")

sona<-function(x){
f<-x
f[!is.na(f)]<-0
f[is.na(f)]<-1
s<-sum(f==1)

f1<-0;
if(f[1]==1)v1=1 else v1=0
for(i in 1:length(f))
if(i<(s+1))
{if(f[1]==1 & f[i]==1 & f[i+1]==1 & f1==0) v1=v1+1 else f1=1}

f2<-0;
if(f[length(f)]==1)v2=length(f) else v2=length(f)+1
for(i in length(f):1)
if(i>(length(f)-s-1)) 
{if(f[length(f)]==1 & f[i]==1 & f[i-1]==1 & f2==0) v2=v2-1 else f2=1}

so.na<-c(v1,length(f)-v2+1,v1,v2)
so.na
}

so.na1<-sona(x)

kt<-0
for(i in 1:length(x))if(i>so.na1[3] && i<so.na1[4] & is.na(x[i]))kt=kt+1
if(kt>0)stop("Series has NA value!")

x1<-x[(so.na1[3]+1):(so.na1[4]-1)]

#Xac dinh khoang truot
sodu=L-(L%/%2)*2
if (sodu==1) a="so le"
else
if (sodu==0) a="so chan"
else
stop("L khong phu hop")

tog<-function(x,L1,L2){
S=0
for(i in L1:L2) S=S+x[i]
S
}

n=length(x1)
F<-c(1:n)

#Tinh toan
#L le
if(a=="so le"){
for(t in 1:n){
if (t<(1+(L-1)/2)) F[t]=NA
else
if (t>(n-(L-1)/2)) F[t]=NA
else
F[t]<-(1/L)*tog(x1,t-(L-1)/2,t+(L-1)/2)
}
}


#L chan

if(a=="so chan"){
F1<-c(1:n)
for(t in 1:n){
if(t<(L/2)) F1[t]=NA
else
if(t>(n-(L/2))) F1[t]=NA
else
F1[t]<-(1/L)*tog(x1,t-(L/2)+1,t+(L/2))
}

F2<-c(1:n)
for(t in 1:n){
if(t<((L/2)-1)) F2[t]=NA
else
if(t>(n-(L/2)-1)) F2[t]=NA
else
F2[t]<-(1/L)*tog(x1,t-(L/2)+2,t+(L/2)+1)
}
F=(F1+F2)/2
}

KQ<-F
KQ
}
