"Kurt" <-
function(x)
{
n <- length (x[!(is.na(x))])
sd <- sqrt(var(x,na.rm=TRUE))
m <- mean(x,na.rm=TRUE)
ku <-(((n*(n+1))/((n-1)*(n-2)*(n-3)))*(sum(((x-m)/sd)^4, na.rm=TRUE)))-((3*(n-1)^2)/((n-2)*(n-3)))
se.ku <- sqrt(24/n)
t.ku <- ku/se.ku
p.ku <- 1-pnorm(abs(t.ku))#Note - evaluated against z instead of t
matK <- cbind(ku,se.ku,t.ku,p.ku)
return(matK)
}

