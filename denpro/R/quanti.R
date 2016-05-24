quanti<-function(values,lkm,base){
#Quantises a vecor of values
#
#values is len-vector
#lkm is positive integer
#base>0
#
#returns len-vector
#
ma<-max(values)
askel<-ma/(lkm-1)
len<-length(values)
ans<-matrix(0,len,1)
for (i in 1:len){
  inv<-base^(values[i]*log(ma+1,base)/ma)-1
  ind<-round(inv/askel)+1
  diskr<-ma*seq(0,lkm-1,1)/(lkm-1)
  disinv<-diskr[ind]
  ans[i]<-ma*log(disinv+1,base)/log(ma+1,base)
}
return(ans)
}
