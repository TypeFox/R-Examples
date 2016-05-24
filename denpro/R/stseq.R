stseq<-function(N,lnum,
refe=NULL,func=NULL,dendat=NULL,
h=NULL,Q=NULL,kernel="epane",weights=NULL,
sig=rep(1,length(N)),support=NULL,theta=NULL,
M=NULL,p=NULL,mul=3,
t=rep(1,length(N)),marginal="normal",r=0,
mu=NULL,xi=NULL,Omega=NULL,alpha=NULL,df=NULL,g=1,
base=10
)
{
#lnum<-length(lseq)
level<-matrix(0,lnum,1)
volume<-matrix(0,lnum,1)
if (!is.null(dendat)) 
  pcf<-pcf.kern(dendat,h,N,kernel=kernel,weights=weights)
else
  pcf<-pcf.func(func,N,   #eval.func.dD
  sig=sig,support=support,theta=theta,
  M=M,p=p,mul=mul,
  t=t,marginal=marginal,r=r, 
  mu=mu,xi=xi,Omega=Omega,alpha=alpha,df=df,g=g)

maksi<-max(pcf$value)
l1<-maksi/(lnum+1) 
lmax<-maksi*lnum/(lnum+1)
level<-hgrid(l1,lmax,lnum,base=base)
level<-level[length(level):1]

for (i in 1:lnum){   
      #lev<-maksi*i/(lnum+1) 
      #level[i]<-lev 
      lev<-level[i]
      if (is.null(refe)) refe<-locofmax(pcf)
      st<-leafsfirst(pcf,lev=lev,refe=refe)
      volume[i]<-max(st$volume)
      if (i==1){
           if (lnum==1){ 
               istseq<-st
           }
           else{
               stseq<-list(st)
           }
      }
      else{
          stseq<-c(stseq,list(st))
      }
}
return(list(shtseq=stseq,level=level,volume=volume,pcf=pcf))
}



