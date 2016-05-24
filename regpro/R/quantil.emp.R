quantil.emp<-function(y,p)
{
n<-length(y)
or<-order(y)
m<-ceiling(p*n)
quant<-y[or[m]]

return(quant)
}


