OOBIndex <-
function(mySample)
{ # return the oob sample,  return vector length=n, if i == 1 , OOB else nonOOB
    n_data<-length(mySample)
    myOOB<-rep(1,n_data)
    myOOB[mySample]<-0
    myOOB<-which(myOOB==1)

   output<-myOOB
}
