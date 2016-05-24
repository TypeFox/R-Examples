ARnames2 <-
function(p,type)
{
if(p==1) tem1 <- "AR"
if(p>1) tem1 <- c("AR",paste("phi",1:(p-1),sep=""))
if(type=="const") tem2 <- c(tem1,"const")
if(type=="const+trend") tem2 <- c(tem1,"const","trend")
return(tem2)
}
