selectand <-
function(df,type=1,step=100,ncol=0,from=0,to=1,col=2){
# Selecting utility
 if (ncol>0){
  if (step<1) step<-100
  n<-dim(df)[1]
  m<-dim(df)[2]
  if ((type==1) | (type==2) | (type==4)){
   xmin<-(-pi)
   xmax<-pi
  } else
  if (type==3){
   xmin<-0
   xmax<-4*pi
  }
  from<-normalize(c(min(df[,ncol]),from,max(df[,ncol])))[2]
  to<-normalize(c(min(df[,ncol]),to,max(df[,ncol])))[2]
  for (i in 1:m) df[,i]<-normalize(df[,i])
  dfm<-numarray(df)
  m<-dim(dfm)[2]
  coorx<-0:step
  for (p in 0:step) coorx[p+1]<-(xmin+(xmax-xmin)/step*p)
  coory<-0:step
  for (i in 1:n) if (df[i,ncol]>=from & df[i,ncol]<=to) {
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
   lines(coorx,coory,col=col)
  }
 }
}

