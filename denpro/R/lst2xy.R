lst2xy<-function(lst,type="radius",gnum=1000)
{
# gives the x and y vectors of a volume transform

if (type=="radius") pv<-plotvolu(lst,data=T,toplot=F)
else{
   lst2<-lst
   lst2$volume<-lst$proba
   pv<-plotvolu(lst2,data=T,toplot=F)
}

lenni<-length(pv$xcoor)/2
xs<-t(matrix(pv$xcoor,2,lenni))
ys<-matrix(0,lenni,1)
for (i in 1:lenni) ys[i]<-pv$ycoor[2*i]

or<-order(ys)
xs<-xs[or,]
ys<-ys[or]

xlow<-min(xs)
xhig<-max(xs)
xstep<-(xhig-xlow)/gnum
x<-seq(xlow,xhig,xstep)
y<-matrix(0,length(x),1)

x[1]<-xs[1,1]
x[length(x)]<-xs[1,2]
i<-2
while (i <= lenni){
  lowind<-round(length(x)*(xs[i,1]-xlow)/(xhig-xlow))
  higind<-round(length(x)*(xs[i,2]-xlow)/(xhig-xlow))
  y[higind:lowind]<-ys[i]
  i<-i+1
}

return(list(x=x,y=y))
}



