cv.filter<-function(data, cutoff=0.05)
{
cv<-apply(data,1,function(x) sd(x)/mean(x))
co<-cutoff
plot(sort(cv), ylab="Coefficient of variation")
abline(h=co,col="red")
data<-data[cv>cutoff,]
return(data)
}