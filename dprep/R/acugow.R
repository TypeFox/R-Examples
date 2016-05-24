acugow <-
function (x,data,vnom=NULL) 
{p=dim(data)[2]
 n=dim(data)[1]
 if (sum(is.na(data))> 0) 
     stop("This dataset has missing values, impute them before running this function.\n",call.=FALSE)
 f=p
 if(length(vnom)>0){fc=(1:f)[-vnom]}
 else {fc=(1:f)}
#print(fc)
dh=rep(1,f)
data1=rbind(x,data)
for (j in fc)
{dh[j]=diff(range(data1[,j],na.rm=TRUE))
    #dh[j]=diff(range(data[,j],na.rm=TRUE))
}
#print(dh)
distall=matrix(0,n,f)
data1=data.matrix(data1)
#print(data1[1,fc])
tempo=sweep(abs(data1[1,fc]-t(data1[-1,fc])),1,dh[fc],"/")
if(length(vnom)>0)
    {for(i in vnom)
    {
    for(j in 2:(n+1))
    {
        if(data1[1,i]!=data1[j,i]){distall[(j-1),i]=1}
    }
}
}
distall[,fc]=t(tempo)
distancia=drop(rowSums(distall)/f)
return(list(matdist=distall,dist=distancia))
}
