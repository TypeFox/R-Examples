massat<-function(rec){
#Calculates a vector of masses of a set of rectangles
#
#rec is k*(2*d)-matrix, represents k rectangles in d-space
#Returns a k-vector
#
#if (dim(t(rec))[1]==1) k<-1 else k<-length(rec[,1])  #rows of rec
if (dim(t(rec))[1]==1){
 d<-length(rec)/2
 vol<-1
 j<-1
   while ((j<=d) && (vol>0)){
     if (rec[2*j]<=rec[2*j-1]) vol<-0
     else vol<-vol*(rec[2*j]-rec[2*j-1])
     j<-j+1
   }
  tulos<-vol
}
else{
 k<-length(rec[,1])
 d<-length(rec[1,])/2
 tulos<-matrix(0,k,1)
 for (i in 1:k){
   vol<-1
   j<-1
   while ((j<=d) && (vol>0)){
     if (rec[i,2*j]<=rec[i,2*j-1]) vol<-0
     else vol<-vol*(rec[i,2*j]-rec[i,2*j-1])
     j<-j+1
   }
   tulos[i]<-vol
 }
}
return(tulos)
}

