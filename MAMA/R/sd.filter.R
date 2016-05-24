sd.filter<-function(data, cutoff=0.5)
{
sd<-apply(data,1,sd)
plot(sort(sd), ylab="Standard deviation")
abline(h=0.5,col="red")
data<-data[sd>cutoff,]
return(data)
}

