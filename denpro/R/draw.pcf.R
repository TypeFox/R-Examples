draw.pcf<-function(pcf,pnum=rep(32,length(pcf$N)),corona=5,dens=FALSE,minval=0,
drawkern=TRUE)
{
#Makes data for drawing a perspective plot.
#pnum on kuvaajan hilan pisteiden lkm
#corona makes corona around the support (useful for densities)

d<-length(pcf$N)

if (d==2){

#col=matrix("black",length(pcf$value),1)

step<-matrix(0,d,1)
for (i in 1:d){
   step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]
}

if ((drawkern)&&(!is.null(pcf$index))){
   return(draw.kern(pcf$value,pcf$index,pcf$N,pcf$support,minval=minval))
}

else{
     pit<-matrix(0,d,1)
     for (i in 1:d){
         pit[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pnum[i]
     }
     alkux<-pcf$support[1]+pit[1]/2-corona*pit[1]
     alkuy<-pcf$support[3]+pit[2]/2-corona*pit[2]
     loppux<-pcf$support[2]-pit[1]/2+corona*pit[1]
     loppuy<-pcf$support[4]-pit[2]/2+corona*pit[2]

     pnum2<-pnum+2*corona
     x<-alkux+c(0:(pnum2[1]-1))*pit[1]
     y<-alkuy+c(0:(pnum2[2]-1))*pit[2]

     reclkm<-length(pcf$value)
     xdim<-length(x)
     ydim<-length(y)
     arvot<-matrix(minval,xdim,ydim)

     l<-1
     while (l<=reclkm){
         begx<-pcf$support[1]+step[1]*(pcf$down[l,1])   #recs[l,1]
         endx<-pcf$support[1]+step[1]*(pcf$high[l,1])     #recs[l,2]
         begy<-pcf$support[3]+step[2]*(pcf$down[l,2])   #recs[l,3]
         endy<-pcf$support[3]+step[2]*(pcf$high[l,2])     #recs[l,4]

         begxind<-round(pnum2[1]*(begx-alkux)/(loppux-alkux))
         endxind<-round(pnum2[1]*(endx-alkux)/(loppux-alkux))
         begyind<-round(pnum2[2]*(begy-alkuy)/(loppuy-alkuy))
         endyind<-round(pnum2[2]*(endy-alkuy)/(loppuy-alkuy))

         arvot[begxind:endxind,begyind:endyind]<-pcf$value[l]
         #col[(begxind+(ydim-1)*begxind):(endxind+(ydim-1)*endyind)]<-ts[l]

         l<-l+1
     }

     #return(apu)
     return(list(x=x,y=y,z=arvot))
}  # else

}  # if (d==2)

else{  #d==1

if (dens){ 
  N<-pcf$N+2 
  alku<-1
}
else{
  N<-pcf$N
  alku<-0
}

x<-matrix(0,N,1)
y<-matrix(0,N,1)
step<-(pcf$support[2]-pcf$support[1])/pcf$N
minim<-pcf$support[1]
maxim<-pcf$support[2]
indenum<-length(pcf$value)

i<-1
while (i<=indenum){

   inde<-pcf$high[i]
   point<-minim+step*inde-step/2
    
   y[alku+inde]<-pcf$value[i]
   x[alku+inde]<-point

   i<-i+1
}

if (dens){
  x[1]<-minim-step/2
  x[N]<-maxim+step/2
}

# remove zeros
if (dens){
  numposi<-0
  for (i in 1:length(y)){
    if (y[i]>0){
       numposi<-numposi+1
       y[numposi]<-y[i]
       x[numposi]<-x[i]
    }
  }
  x<-x[1:numposi]
  y<-y[1:numposi]
}

or<-order(x)
xor<-x[or]
yor<-y[or]

return(list(x=xor,y=yor))

}  # else d=1

}











