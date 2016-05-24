SGoF <-
function(u,alpha=0.05,gamma=0.05){

sgof<-function(u,alpha=0.05,gamma=0.05){

v=as.numeric(u<=gamma)
n=length(v)


SGoF=max(min(floor(n*(mean(v)-gamma)-n*sqrt(mean(v)*(1-mean(v))/n)*qnorm(1-alpha)+1),sum(as.integer(n*ecdf(u)(u))<=n*(mean(v)-gamma)-n*sqrt(mean(v)*(1-mean(v))/n)*qnorm(1-alpha)+1)),0)

su<-sort(u)
jj<-which(u==1)
if(length(jj)!=0) pi0<-1 else pi0<-min((-1/n)*sum(log(1-u)),1)

if(SGoF==0){FDR_S<-0}else{FDR_S<-round((pi0*su[SGoF])/(ecdf(u)(su[SGoF])),4)}


Nu1=pmax(n*(ecdf(su)(su)-su)-sqrt(n*ecdf(su)(su)*(1-ecdf(su)(su)))*qnorm(1-su)+1,0)


Nu2<-sapply(1:n,function(i) Nu2<-max(which(as.integer(n*ecdf(su)(su))<=Nu1[i]),0))

Nu<-pmin(Nu1,Nu2)



Nmax=max(Nu,na.rm=T)


a.p=rep(1,n)  
for (i in 1:Nmax) {a.p[i]=min(su[as.integer(n*ecdf(su)(su[i]))<=Nu],na.rm=T)}

return(c(list(Rejections=SGoF,FDR=min(FDR_S,1),Adjusted.pvalues=a.p)))
}






if(missing(u)){stop("data argument is required")}
u<-as.vector(u)
res<-sgof(u,alpha,gamma)
res$data<-sort(u)
res$alpha<-alpha
res$gamma<-gamma
res$call<-match.call()
class(res)<-"SGoF"
return(res)
}
