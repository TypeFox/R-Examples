`bray.curtis` <-
function(x,y) 
{
bray<-1-sum(abs(x-y))/sum(x+y)
return(bray)
}

