cenone<-function(rec){
#Calculates the 1st moment of a rectangle.
#
#rec is (2*d)-vector, represents rectangle in d-space
#Returns a d-vector.
#
d<-length(rec)/2
res<-matrix(0,d,1)
for (j in 1:d){
  apurec<-rec      #apurec such that is volume is equal to
  apurec[2*j-1]<-0 #volume of d-1 dimensional rectangle where
  apurec[2*j]<-1   #we have removed j:th dimension
  vajmas<-massone(apurec) 
  res[j]<-vajmas*(rec[2*j]^2-rec[2*j-1]^2)/2  
}
return(res) 
}

