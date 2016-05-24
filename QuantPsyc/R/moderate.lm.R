"moderate.lm" <-
function (x,z,y,data,mc=FALSE)
{
attach(data)
if(!mc)
{mcx <- x - mean(x, na.rm=TRUE)
mcz <- z - mean(z, na.rm=TRUE)}
else {
mcx <- x
mcz <- z }

lm1 <- lm(y ~ mcx*mcz, na.action=na.omit)
detach(data)
return(lm1)
}

