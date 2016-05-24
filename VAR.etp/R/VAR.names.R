VAR.names <-
function(x,p,type="const")
{
tem1 <- colnames(x)
varnames <- character()
for(i in 1:p)
{
tem2 <- paste(tem1,rep(-i,ncol(x)),sep="(")
tem3 <- paste(tem2,")",sep="")
varnames <- c(varnames,tem3)
}
if(type=="const") varnames <- c(varnames,"const")
if(type=="const+trend") varnames <- c(varnames,"const","trend")
return(varnames)
}
