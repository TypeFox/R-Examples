selecto <-
function(df,ncol=0,from=0,to=1,col=2){
# Selecting utility
 if (ncol>0){
  n<-dim(df)[1]
  m<-dim(df)[2]
  from<-normalize(c(min(df[,ncol]),from,max(df[,ncol])))[2]
  to<-normalize(c(min(df[,ncol]),to,max(df[,ncol])))[2]
  for (i in 1:m) df[,i]<-normalize(df[,i])
  numattr<-nnumattr(df)
  coory<-1:numattr
  for (i in 1:n) if (df[i,ncol]>=from & df[i,ncol]<=to) {
   naxis<-numattr
   for (a in 1:m){
    if (is.numeric(df[,a])){
     coory[numattr-naxis+1]<-df[i,a]
     naxis<-naxis-1
    }
   }
   lines(1:numattr,coory,col=col)
  }
 }
}

