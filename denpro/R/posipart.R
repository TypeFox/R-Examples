posipart<-function(pcf)
{
pcf$value<-pmax(pcf$value,0)
return(pcf)
}

