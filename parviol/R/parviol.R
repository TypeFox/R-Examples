parviol <-
function(df,violinplot=TRUE,main=NULL,sub=NULL){
 # Setting plot
 n<-dim(df)[1]
 m<-dim(df)[2]
 plot.new()
 numattr<-nnumattr(df)
 plot.window(c(0.5,numattr+0.5),c(-0.2,1.2))
 title(main=main,sub=sub)
 # Asix
 naxis<-numattr
 for (i in 1:m){
  if (is.numeric(df[,i])){
   lines(c(numattr-naxis+1,numattr-naxis+1),c(0,1))
   text(numattr-naxis+1,1.2,colnames(df)[i],cex=0.7)
   text(numattr-naxis+1,1.05,max(df[,i]),cex=0.6)
   text(numattr-naxis+1,-0.05,min(df[,i]),cex=0.6)
   naxis<-naxis-1
  }
 }
 # Normalization
 for (i in 1:m) df[,i]<-normalize(df[,i])
 # Violin plots
 naxis<-numattr
 if (violinplot) for (i in 1:m){
  if (is.numeric(df[,i])){
   vioplot(df[,i],add=TRUE,at=numattr-naxis+1,col=8,border=8,rectCol=rgb(0.3,0.3,0.3))
   naxis<-naxis-1
  }
 }
 # Parallel coordinates
 coory<-1:numattr
 for (i in 1:n){
  naxis<-numattr
  for (a in 1:m){
   if (is.numeric(df[,a])){
    coory[numattr-naxis+1]<-df[i,a]
    naxis<-naxis-1
   }
  }
  lines(1:numattr,coory)
 }
}

