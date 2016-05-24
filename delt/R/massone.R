massone<-function(rec){
#Calculates the mass of a rectangle.
#
#rec is (2*d)-vector, represents rectangle in d-space
#Returns a real number >0.
#
d<-length(rec)/2
vol<-1
for (j in 1:d){
  vol<-vol*(rec[2*j]-rec[2*j-1])
}
return(vol) 
}

