disc.1r <-
function (data,convar,binsize=15,out=c("symb","num")) 
{#data=as.matrix(data)
if(sum(is.na(data))> 0) 
stop("This dataset has missing values, impute them before running this function.\n",call.=FALSE)
    p=dim(data)[2]
n=dim(data)[1]
data1=data
data1[,p]=as.numeric(factor(data1[,p]))
#print(data1)
       for (i in convar)
{    a=data1[,c(i,p)]
     
#b=a[order(a[,1]),]
data1[,i]=unor(a,binsize,out)
data[,1]=data[,1]}
as.data.frame(cbind(data1[,-p],data[,p]))
}
