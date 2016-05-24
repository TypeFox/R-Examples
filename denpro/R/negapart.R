negapart<-function(pcf)
{
pcf$value<--pmin(pcf$value,0)
return(pcf)
}
