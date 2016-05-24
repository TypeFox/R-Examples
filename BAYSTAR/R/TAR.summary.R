TAR.summary <-function(x,lagp1,lagp2,constant=1)
{
n<-ncol(x)
temp<-matrix(NA,n,5)
for(i in 1:n){
temp[i,1]<-mean(x[,i])
temp[i,2]<-quantile(x[,i],0.5)  ## Median of estimates
temp[i,3]<-sd(x[,i])
temp[i,4:5]<-quantile(x[,i],c(0.025,0.975))    ## 95% Bayes interval of estimates

colnames(temp)<-c("mean","median","s.d.","lower","upper")
if(constant==1)
{
rownames(temp)<-c(paste("phi1",c(0,lagp1),sep="."),paste("phi2",c(0,lagp2),sep="."),"sigma^2 1","simga^2 2","r","mean1","mean2")
}
else
{
rownames(temp)<-c(paste("phi1",lagp1,sep="."),paste("phi2",lagp2,sep="."),"sigma^2 1","simga^2 2","r")
}
}
return(temp)
}

