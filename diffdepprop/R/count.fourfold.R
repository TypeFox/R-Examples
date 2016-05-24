count.fourfold=function(data,col.test1,col.test2)
{
a=sum(data[,col.test1]==1 & data[,col.test2]==1)
b=sum(data[,col.test1]==1 & data[,col.test2]==0)
c=sum(data[,col.test1]==0 & data[,col.test2]==1)
d=sum(data[,col.test1]==0 & data[,col.test2]==0)
n=a+b+c+d
cint=c(a,b,c,d, n)
return(list(paste("The elements of the fourfold table are: a=",a,", b=",b,", c=",c,", d=",d, ". The total sample size is: n=",n),cint))

}
