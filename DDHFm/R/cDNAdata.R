"cDNAdata" <-
function(data.vect,cdnalength,datasize,ng,nrep){


##Reads the data

data1<-matrix(data.vect,ng,nrep)


#samples from the full data

data1<-data1[1:cdnalength,]
s1<-seq(1,cdnalength,1)
sample1<-sample(s1,datasize)
datamat<-data1[sample1,]

means<-matrix(0,datasize,1)
sds<-matrix(0,datasize,1)

for(i in 1:datasize){
means[i,]<-mean(datamat[i,])
sds[i,]<-sqrt(var(datamat[i,]))
}

means<-c(means)

smeans<-sort.list(means)
sorteddata<-datamat[smeans,]

sorteddata.y<-t(sorteddata)
y<-c(sorteddata.y)

return(y)
}

