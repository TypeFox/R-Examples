ML.k<-function(f,x,res=0.000001){
N<-sum(f)
m<-sum(x*(f/N))
Sfx<-sum(f*x)
Sfx2<-sum(f*x^2)
s2<-(Sfx2-(Sfx^2/N))/(N-1)
khat1<-m^2/(s2-m)
rng<-s2/10

Az<-matrix(ncol=1,nrow=length(f))
Az[1]<-N-f[1]
for(i in 2:length(f)){
Az[i]<-Az[i-1]-f[i]
}

khats<-seq(khat1-rng,khat1+rng,by=res)
z<-matrix(ncol=1,nrow=length(khats))
for(i in 1:length(khats)){
S1<-sum(Az/(khats[i]+x))
S2a <- 1+(m/khats[i])
S2<-ifelse(S2a >= 0, N*log(S2a), 0)
z[i]<-abs(S1-S2)
}

result<-list()
result$k<-khats[z==min(z)]
result$m<-m
result
}
