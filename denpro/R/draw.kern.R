draw.kern<-function(value,index,N,support,minval=0,dendat=NULL,h=NULL)
{

d<-length(N)

if (d==2){

x<-matrix(0,N[1]+2,1)
y<-matrix(0,N[2]+2,1)
z<-matrix(minval,N[1]+2,N[2]+2)
#col<-matrix("black",dim(z)[1]*dim(z)[2],1)

minim<-matrix(0,d,1)  #minim is d-vector of minimums
maxim<-matrix(0,d,1)
for (i in 1:d){
  minim[i]<-support[2*i-1]
  maxim[i]<-support[2*i]
}
delta<-(maxim-minim)/(N+1)

indenum<-dim(index)[1]

i<-1
while (i<=indenum){
   inde<-index[i,]
   z[1+inde[1],1+inde[2]]<-value[i]
   #col[1+inde[1]+dim(z)[1]*inde[2]]<-ts[i]
   i<-i+1
}

i<-1
while (i<=N[1]){
   x[1+i]<-support[1]+delta[1]*i
   i<-i+1
}

i<-1
while (i<=N[2]){
   y[1+i]<-support[3]+delta[2]*i
   i<-i+1
}

x[1]<-support[1]
x[N[1]+2]<-support[2]
y[1]<-support[3]
y[N[2]+2]<-support[4]

return(list(x=x,y=y,z=z)) #col=col[length(col):1]))

}

else{    #d=1

x<-matrix(0,N+2,1)
y<-matrix(0,N+2,1)

minim<-min(dendat)
maxim<-max(dendat)
hmax<-max(h)
delta<-(maxim-minim+2*hmax)/(N+1)

indenum<-dim(index)[1]

i<-1
while (i<=indenum){

   inde<-index[i]
   point<-minim-hmax+delta*inde
    
   y[1+inde]<-value[i]
   x[1+inde]<-point

   i<-i+1
}
x[1]<-minim-hmax
x[N+2]<-minim-hmax+delta*N+delta

return(list(x=x,y=y))

}

}











