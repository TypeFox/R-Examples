depths1<-function(m,j){
if(m < j)depths1<-0
else{
if(j==1)depths1<-m
if(j==2)depths1<-(m*(m-1))/2
if(j==3)depths1<-(m*(m-1)*(m-2))/6
}
depths1
}
