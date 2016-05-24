BH <-
function(u,alpha=0.05){

bh<-function(u,alpha=0.05){

n=length(u)
r=rank(u,ties.method="max")


bh=max(c(r[u<=(r/n)*alpha],0),na.rm = T) 

su<-sort(u)
jj<-which(u==1)
if(length(jj)!=0) pi0<-1 else pi0<-min((-1/n)*sum(log(1-u)),1)

if(bh==0){FDR_BH<-0}else{FDR_BH<-round((pi0*su[bh])/(ecdf(u)(su[bh])),4)}



ad.p=numeric(n)
ad.p[n]<-sort(u)[n]
for(i in (n-1):1){
ad.p[i]<-min(sort(u)[i]*(n/i),ad.p[i+1])
}


return(c(list(Rejections=bh,FDR=min(FDR_BH,1),Adjusted.pvalues=sort(ad.p))))
}

if(missing(u)){stop("data argument is required")}




u<-as.vector(u)
res<-bh(u,alpha)
res$data<-sort(u)
res$alpha<-alpha
res$call<-match.call()
class(res)<-"BH"
return(res)
}
