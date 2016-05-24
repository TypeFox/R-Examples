`optim.fun.lm` <-
function(x, data, param)
{
l1 <- x[1]
l2 <- x[2]
l3 <- x[3]
l4 <- x[4]
data <- sort(data)
check1 <- gl.check.lambda.alt(l1, l2, l3, l4, param)
check2 <- qgl(0, l1, l2, l3, l4, param) <= min(data)
check3 <- qgl(1, l1, l2, l3, l4, param) >= max(data)
overall<-sum(c(check1,check2,check3),na.rm=T)

if(overall==3){
response<-sum((Lmoments(data)-fun.lm.theo.gld(l1,l2,l3,l4,param))^2)
}
if(overall!=3){
response <- NA
}
return(response)
}

