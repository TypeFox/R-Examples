print.maxstat<-function(x,...)
{
 printCoefmat(t(x), dig.tst=3, tst.ind=1:4, P.values = TRUE, na.print="-", ...)
}