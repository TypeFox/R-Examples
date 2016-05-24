point.eval<-function(tr,x)
{
# tr is an evaluation tree

d<-length(tr$support)/2
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(tr$support[2*i]-tr$support[2*i-1])/tr$N[i]

# if x is not in the support, then ans=0
insupport<-1
for (i in 1:d){
    if ((x[i]>=tr$support[2*i]) || (x[i]<=tr$support[2*i-1])){
       ans<-0
       insupport<-0
    }
}
if (insupport==1){
  node<-1
  while (tr$left[node]>0){
      dir<-tr$direc[node]
      spl<-tr$split[node]
      realspl<-tr$support[2*dir-1]+spl*step[dir]
      if (x[dir]>realspl) node<-tr$right[node]
      else node<-tr$left[node]
  }
  #loc<-tr$infopointer[node]
  #ans<-tr$value[loc]
  ans<-tr$mean[node]
}

return(ans)
}
