"qqnormy" <-
function (y) 
{

v<-qqnorm(y,plot.it=FALSE)
return(v$y[sort.list(v$x)])

}

