AR.order <-
function(x,pmax)
{
n <- nrow(x)
tem <- numeric(0)
for( i in 1:pmax)
{
e <- LSM(x,i)$resid
aic <- log( sum(e^2)/n ) + (2/n) * (i+1)
bic <- log( sum(e^2)/n )+ (log(n)/n) * (i+1)
hq <- log(sum(e^2)/n ) + (2*log(log(n))/n) * (i+1)

tem <- rbind(tem,cbind(i,aic,bic,hq))
}
paic <- tem[tem[,2] == min(tem[,2]),1]
pbic <- tem[tem[,3] == min(tem[,3]),1]
phq <- tem[tem[,4] == min(tem[,4]),1]

return(list(ARorder=cbind(paic,pbic,phq),Criteria=tem[,2:4]))
}
