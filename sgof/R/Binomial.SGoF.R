Binomial.SGoF <-
function(u,alpha=0.05,gamma=0.05){

binomial.sgof<-function(u,alpha=0.05,gamma=0.05){

v=as.numeric(u<=gamma)
n=length(v)

b=seq(0,n)



Binomial.SGoF=max(min(n*ecdf(u)(gamma)-b[min(which((1-pbinom(b-1,n,gamma))<=alpha))]+1,sum(as.integer(n*ecdf(u)(u))<=(n*ecdf(u)(gamma)-b[min(which((1-pbinom(b-1,n,gamma))<=alpha))]+1))),0)


su<-sort(u)
jj<-which(u==1)
if(length(jj)!=0) pi0<-1 else pi0<-min((-1/n)*sum(log(1-u)),1)

if(Binomial.SGoF==0){FDR_BS<-0}else{FDR_BS<-round((pi0*su[Binomial.SGoF])/(ecdf(u)(su[Binomial.SGoF])),4)}

su=sort(u)

cuantiles<-sapply(1:n,function(i) cuantiles<-b[min(which((1-pbinom(b-1,n,su[i]))<=su[i]))])


N<-pmax(n*ecdf(su)(su)-cuantiles+1,0)




N1<-sapply(1:n,function(i) N1<-max(which(as.integer(n*ecdf(su)(su))<=N[i]),0))



N2<-pmin(N,N1)



Nmax=max(N2[N2!=n])
auu=rep(1,n)  

for (i in 1:Nmax) {auu[i]=min(su[as.integer(n*ecdf(su)(su[i]))<=N2])}


return(c(list(Rejections=Binomial.SGoF,FDR=min(FDR_BS,1),Adjusted.pvalues=auu)))
}




if(missing(u)){stop("data argument is required")}
u<-as.vector(u)
res<-binomial.sgof(u,alpha,gamma)
res$data<-sort(u)
res$alpha<-alpha
res$gamma<-gamma
res$call<-match.call()
class(res)<-"Binomial.SGoF"
return(res)
}
