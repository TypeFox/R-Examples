"is.inf" <-
function(x)
{
result<-(x==Inf)+(x==-Inf)
result<-as.logical(result)

return(result)
}

