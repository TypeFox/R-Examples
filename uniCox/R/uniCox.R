
uniCox=function(x,y,status,lamlist=NULL,nlam=20,del.thres=.01, max.iter=5){

# x is n by p
 this.call <- match.call()
if(nrow(x)!=length(y)){
           stop("Error in data input")
         }
if(nrow(x)!=length(status)){
           stop("Error in data input")
         }

if((sum(status==0)+sum(status==1))!=length(status)){
          stop("Error in data input- status")
         }
 

        


mx=colMeans(x)
x=scale(x,center=mx,scale=F)
vx=1/coxvar(t(x), y, status)
s0=quantile(sqrt(vx),.5)
xs=scale(x,center=F,scale=sqrt(vx)+s0)

u0=coxscor(t(xs),y,status)$scor
inffull=NULL
if(is.null(lamlist)){lamlist=seq(0,max(abs(u0)),length=nlam)}
nlam=length(lamlist)

p=ncol(xs)
beta2=matrix(0,nrow=p,ncol=length(lamlist))
beta0=rep(0,p)
for(k in 1:length(lamlist)){
cat(c("lambda value ",k,"out of ",length(lamlist)),fill=T)
lam=lamlist[k]
 om=abs(u0)>lam
 beta=rep(0,sum(om))
 if(k>1){beta=beta2[om,k-1,drop=F]}
niter=0
u=rep(1,p)
go=T
  while(niter< max.iter & go &sum(om)){
   niter=niter+1
   offset=t(scale(xs[,om,drop=F],center=F,scale=1/beta))
    u=coxscor2(t(xs[,om,drop=F]),y,status,offset=offset)$scor-lam*sign(beta)
    v=coxvar2(t(xs[,om,drop=F]), y, status,offset=offset)
    obeta=beta
    beta=obeta+(u/v)
    del=mean(abs(beta-obeta),na.rm=TRUE)
   if(del< del.thres){go=F}
}
if(lam==0){betafull=beta;inffull=v}
beta2[om,k]=beta
}

nfeatures=colSums(beta2!=0, na.rm=T)
nm=sum(is.na(beta2))
if(nm>0){cat(c(nm," betas missing"),fill=T)}
junk=list(lamlist=lamlist,beta=beta2,nfeatures=nfeatures,mx=mx,vx=vx,s0=s0, call=this.call)
class(junk)="uniCoxFit"
return(junk)
}


