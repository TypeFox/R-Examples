drawgene<-function(values,recs,plkm=60,ep1=0.5){
#Makes data for drawing a perspective plot.

#plkm on kuvaajan hilan pisteiden lkm
#ep1 makes 0-corona around the support (useful for densities)

#koe<-drawgene(values,recs,plkm=30)
#persp(koe$x,koe$y,koe$z,phi=30,theta=60) 

alkux<-min(recs[,1])-ep1
alkuy<-min(recs[,3])-ep1
loppux<-max(recs[,2])+ep1
loppuy<-max(recs[,4])+ep1
pitx<-(loppux-alkux)/plkm
pity<-(loppuy-alkuy)/plkm
x<-alkux+c(0:plkm)*pitx
y<-alkuy+c(0:plkm)*pity

reclkm<-length(values)
xdim<-length(x)
ydim<-length(y)
arvot<-matrix(0,xdim,ydim)

l<-1
while (l<=reclkm){
   begx<-recs[l,1]
   endx<-recs[l,2]
   begy<-recs[l,3]
   endy<-recs[l,4]

   begxind<-round(plkm*(begx-alkux)/(loppux-alkux))
   endxind<-round(plkm*(endx-alkux)/(loppux-alkux))
   begyind<-round(plkm*(begy-alkuy)/(loppuy-alkuy))
   endyind<-round(plkm*(endy-alkuy)/(loppuy-alkuy))

   arvot[begxind:endxind,begyind:endyind]<-values[l]

   l<-l+1
}

return(list(x=x,y=y,z=arvot))
#persp(x,y,arvot)
}










