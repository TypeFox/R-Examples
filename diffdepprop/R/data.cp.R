data.cp=function(dat,col.test1, col.test2, cp.test1,cp.test2)
{
cp1=c()
cp2=c()
n=length(dat[,col.test1])
for (i in 1:n){
if (dat[i,col.test1]>=cp.test1){cp1[i]=1}
if (dat[i,col.test1]<cp.test1){cp1[i]=0}
if (dat[i,col.test2]>=cp.test2){cp2[i]=1}
if (dat[i,col.test2]<cp.test2){cp2[i]=0}
}
data_out<-data.frame(cp1,cp2)
return(data_out)
}

