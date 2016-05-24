lokeroi<-function(A,b){
#
#a is a vector 0<=A_1<...<A_m
#b is number >0
#
#we want to find index i in 1,...,m so that A_{i} < b <= A_i  
#huom jos i=m, niin b>A_{m}
#
m<-length(A)
res<-m
i<-m-1
while (i>=1){
  if (b<=A[i]) res<-i
  i<-i-1
}
#
return(res)
}










