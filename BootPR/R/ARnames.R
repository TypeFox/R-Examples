ARnames <-
function(p,type)
{tem1 <- paste("AR",1:p,sep="")
if(type=="const") tem2 <- c(tem1,"const")
if(type=="const+trend") tem2 <- c(tem1,"const","trend")
return(tem2)
}
