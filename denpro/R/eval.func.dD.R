eval.func.dD<-function(func,N,
sig=rep(1,length(N)),support=NULL,theta=NULL,g=1,
M=NULL,p=NULL,mul=3,
t=rep(1,length(N)),marginal="normal",r=0,
mu=NULL,xi=NULL,Omega=NULL,alpha=NULL,df=NULL,a=0.5,b=0.5
)   
# func== "mixt", "epan", "cop1"
{
d<-length(N)
recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

if (is.null(support)){

   if (func=="mixt"){
     support<-matrix(0,2*d,1)
        for (i in 1:d){
           support[2*i-1]<-min(M[,i]-mul*sig[,i])
           support[2*i]<-max(M[,i]+mul*sig[,i])
        }
   }

   if (func=="epan"){
      if (is.null(sig)) sig<-c(1,1)
      support<-matrix(0,2*d,1)
      for (i in 1:d){
          support[2*i-1]<--sig[i]
          support[2*i]<-sig[i]
      }
   }
}

if ((marginal=="unif")) support<-c(0,sig[1],0,sig[2])
# && (is.null(support))) 
#support<-c(-sig[1]/2,sig[1]/2,-sig[2]/2,sig[2]/2)


lowsuppo<-matrix(0,d,1)
for (i in 1:d) lowsuppo[i]<-support[2*i-1]
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

numpositive<-0
for (i in 1:recnum){
    inde<-digit(i-1,N)+1
    #if ((inde[1]==0) && (inde[2]==N[2])) inde<-c(0,0)
    point<-lowsuppo+step*inde-step/2

    if (!is.null(theta)){
         rotmat<-matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
         point<-rotmat%*%point
    }

    if (func=="prod") val<-eva.prod(point,marginal,g)
    if (func=="skewgauss") val<-eva.skewgauss(point,mu,sig,alpha)
    #if (func=="dmsn") val<-dmsn(point,xi,Omega,alpha)
    if (func=="student") val<-eva.student(point,t,marginal,sig,r,df)
    if (func=="gumbel") val<-eva.copula(point,
        type="gumbel",marginal=marginal,sig=sig,r=r,t=t,g=g)
    if (func=="frank") val<-eva.copula(point,
        type="frank",marginal=marginal,sig=sig,t=t,g=g)
    if (func=="plackett") val<-eva.plackett(point,t,marginal,sig)
    if (func=="clayton2") val<-eva.clayton(point,t,marginal,sig,df)
    if (func=="clayton") val<-eva.copula(point,
        type="clayton",marginal=marginal,sig=sig,r=r,t=t,g=g)
    if (func=="cop6") val<-eva.cop6(point,t,marginal,sig)
    if (func=="epan") val<-epan(point)
    if (func=="gauss") val<-eva.copula(point,
        type="gauss",marginal=marginal,sig=sig,r=r,t=t)
    if (func=="normal") val<-eva.gauss(point,t=t,marginal=marginal,sig=sig,r=r)
    if (func=="mixt"){
        val<-0
        mixnum<-length(p)
        for (mi in 1:mixnum){
            evapoint<-(point-M[mi,])/sig[mi,]
            val<-val+p[mi]*evanor(evapoint)/prod(sig[mi,])
        }
    }
   if (func=="hat") val<-eva.hat(point,a=a,b=b)

    if (val>0){
       numpositive<-numpositive+1
       value[numpositive]<-val
       index[numpositive,]<-inde
    }
}

value<-value[1:numpositive]
index<-index[1:numpositive,]
down<-index-1
high<-index

res<-list(
value=value,index=index,
down=down,high=high,  #step=delta,
support=support,N=N)

return(res)
}                              











