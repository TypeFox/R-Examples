#transform dataset to specified m and sd
dtrans<-function(data,m,sd,rnd=F){
#basic checking  
x<-dim(data)[2] #getting number of vars
cat("Number of variables in dataset:",  x,"\n")#number of vars
cat("Number of means specified:",  length(m),"\n")
cat("Number of standard deviations in dataset:",  length(sd),"\n")
  if(length(m) != length(sd)){stop("Number of means/SDs should match.")}
  if(length(m) != dim(data)[2]){stop("Incorrect number of means/SDs")}
for (i in 1:x){  
xrange<-range(data[,1])
if(xrange[1]-xrange[2]==0){stop("Constant value detected")}
}

#transforming
for (i in 1:x){
  data[i]=m[i]+sd[i]*data[i]
 }
if(rnd==T){return(round(data,0))
} else {return(data)}
}
