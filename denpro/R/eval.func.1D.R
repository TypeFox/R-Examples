eval.func.1D<-function(func,N,support=NULL,g=1,std=1,distr=FALSE,
M=NULL,sig=NULL,p=NULL,a=0.5,b=0.5,d=2)
{
if (func=="gauss"){
   norma<-(2*pi)^(-1/2)
   funni<-function(t){ fu<-exp(-t^2/2); return( norma*fu ) }
}
if (func=="polynomial"){
   support<-c(-std,std)
   norma<-(2*(1-1/(g+1)))^(-1)
   funni<-function(t){ fu<-1-abs(t)^g; return( norma*fu ) }
}
if (func=="student"){
   norma<-gamma((g+1)/2)/((g*pi)^(1/2)*gamma(g/2))
   funni<-function(t){ fu<-(1+t^2/g)^(-(g+1)/2); return( norma*fu ) }
   #y<-dt(x,df=g)
}
if (func=="exponential"){
   norma<-1/2
   funni<-function(t){ fu<-exp(-abs(t)); return( norma*fu ) }
}
if (func=="exponential"){
   norma<-1/2
   funni<-function(t){ fu<-exp(-abs(t)); return( norma*fu ) }
}
if (func=="mixt"){
   funni<-function(t){ 
       mixnum<-length(p)
       val<-0
       for (mi in 1:mixnum){
               evapoint<-(t-M[mi])/sig[mi]
               val<-val+p[mi]*evanor(evapoint)/sig[mi]
        }
        return( val ) 
   }
}
if (func=="hat"){
   normavak<-((2*pi)^d*(a^(-d)-b))^(-1)
   norma<-normavak*(2*pi)^((d-1)/2)
   funni<-function(t){  #(t,a,b,d,...){ 
          fu<-a^(1-d)*exp(-a^2*t^2)-b*exp(-t^2/2); return( norma*fu ) }
}

if (is.null(support)) support<-c(-1,1)

value<-matrix(0,N,1)
step<-(support[2]-support[1])/N
lowsuppo<-support[1]

if (!distr){
   for (i in 1:N){
       inde<-i
       t<-lowsuppo+step*inde-step/2
       value[i]<-funni(t/std)/std
   }
}
else{
   inde<-1
   t<-lowsuppo+step*inde-step/2
   value[1]<-step*funni(t/std)/std
   for (i in 2:N){
       inde<-i
       t<-lowsuppo+step*inde-step/2
       value[i]<-value[i-1]+step*funni(t/std)/std
       #funni(t/std,g=g,a=a,b=b,d=d)/std
   }
}

index<-seq(1:N)
len<-length(index)
down<-matrix(0,len,1)
high<-matrix(0,len,1)
down[,1]<-index-1
high[,1]<-index

res<-list(
value=value,
down=down,high=high,
#down=index-1,high=index,  
support=support,N=N)

return(res)
}

                              

