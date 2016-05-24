andrews <-
function(df,type=1,clr=NULL,step=100,ymax=10,main=NULL,sub=NULL){
 # Setting plot
 if (step<1) step<-100
 n<-dim(df)[1]
 m<-dim(df)[2]
 plot.new()
 if ((type==1) | (type==2) | (type==4)){
  xmin<-(-pi)
  xmax<-pi
 } else
 if (type==3){
  xmin<-0
  xmax<-4*pi
 }
 plot.window(c(xmin,xmax),c(-ymax,ymax))
 title(main=main,sub=sub)
 axis(1)
 axis(2)
 box()
 lines(c(xmin,xmax),c(0,0))
 # Normalization
 for (i in 1:m) df[,i]<-normalize(df[,i])
 # Colors
 clx<-rep(1,n)
 if (!is.null(clr)){
  if (!is.numeric(df[,clr])){
   for (i in 1:n){
    for (a in 1:nlevels(df[,clr])) if (levels(df[,clr])[a]==df[i,clr]) clx[i]<-rainbow(nlevels(df[,clr]))[a]
   }
  } else
  {
   for (i in 1:n) clx[i]<-hsv(0,1,df[i,clr])
  }
 }
 # Array
 dfm<-numarray(df)
 m<-dim(dfm)[2]
 # Curves
 coorx<-0:step
 for (p in 0:step) coorx[p+1]<-(xmin+(xmax-xmin)/step*p)
 coory<-0:step
 for (i in 1:n){
  for (p in 0:step){
   coory[p+1]<-0
   tt<-(xmin+(xmax-xmin)/step*p)
   for (a in 1:m){
    if (type==1){
     if (a==1){
      coory[p+1]<-dfm[i,a]/(2^0.5)
     } else
     {
      cnst<-(a-2)%/%2+1
      if ((a-2)%%2==0) coory[p+1]<-coory[p+1]+dfm[i,a]*sin(cnst*tt) else coory[p+1]<-coory[p+1]+dfm[i,a]*cos(cnst*tt)
     }
    } else
    if (type==2){
     cnst<-(a-1)%/%2+1
     if ((a-1)%%2==0) coory[p+1]<-coory[p+1]+dfm[i,a]*sin(cnst*tt) else coory[p+1]<-coory[p+1]+dfm[i,a]*cos(cnst*tt)
    } else
    if (type==3){
     coory[p+1]<-coory[p+1]+dfm[i,a]*cos((a*tt)^0.5)
    } else
    if (type==4){
     if (a==1){
      coory[p+1]<-dfm[i,a]
     } else
     {
      cnst<-(a-2)%/%2+1
      if ((a-2)%%2==0) coory[p+1]<-coory[p+1]+dfm[i,a]*(sin(cnst*tt)+cos(cnst*tt)) else 
                       coory[p+1]<-coory[p+1]+dfm[i,a]*(sin(cnst*tt)-cos(cnst*tt))
     }
     coory[p+1]<-coory[p+1]/(2^0.5)
    }
   }
  }
  lines(coorx,coory,col=clx[i])
 }
}

