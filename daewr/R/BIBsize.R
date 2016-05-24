BIBsize<-function(t,k)
{
  b<-t
  r<-0
  lambda<-0
  check<-0
  while (check==0) {
   while (r==0) {
     #cat("r=",r)
     testr<-(b*k)/t
     #cat("testr=",testr,"b=",b)
     if (testr==floor(testr)) {
       r<-testr
       } else {
       b<-b+1
       }
     }
      #cat("b=",b, "r=",r)
      testl<-(r*(k-1))/(t-1)
      #cat("testl=",testl,"b=",b)
      if (testl==floor(testl)) {
       lambda<-testl
       check=1
       } else {
       r<-0
       b<-b+1
      #cat("b=",b, "r=",r)
        }

    #cat("lambda=",lambda)
    }
  cat("Posible BIB design with b=",b," and r=",r," lambda=",lambda,"\n") 
}

