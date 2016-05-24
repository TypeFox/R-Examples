BY <-
function(u,alpha=0.05){

by<-function(u,alpha=0.05){

n=length(u)
r=rank(u,ties.method="max")
q <- sum(1/(1:n))

by=max(c(r[u<=(r/(n*q))*alpha],0),na.rm = T) 

su<-sort(u)
jj<-which(u==1)
if(length(jj)!=0) pi0<-1 else pi0<-min((-1/n)*sum(log(1-u)),1)

if(by==0){FDR_BY<-0}else{FDR_BY<-round((pi0*su[by])/(ecdf(u)(su[by])),4)}



ad.p=numeric(n)
ad.p[n]<-sort(u)[n]
for(i in (n-1):1){
ad.p[i]<-min(sort(u)[i]*q*n/i,ad.p[i+1])
}


return(c(list(Rejections=by,FDR=min(FDR_BY,1),Adjusted.pvalues=sort(ad.p))))
}

if(missing(u)){stop("data argument is required")}




u<-as.vector(u)
res<-by(u,alpha)
res$data<-sort(u)
res$alpha<-alpha
res$call<-match.call()
class(res)<-"BY"
return(res)
}
