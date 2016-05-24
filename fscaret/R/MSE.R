MSE<-function(vect1,vect2,rows_no){

result=0
pred<-0
obs<-0

for(i in 1:rows_no){
  result<-result +(vect1[i]-vect2[i])^2
  }
  result<-(result/rows_no)
  return(result)
}