a_st <-
function(x, na.rm=F)
{
x<-x-min(x, na.rm=na.rm)
x<-x/max(x, na.rm=na.rm)
x
}
