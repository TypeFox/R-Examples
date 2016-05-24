# aus e1071
countpattern2<-function (x)
{
    intpatt<-apply(as.matrix(x),1,bin2int)
    intpatt<-intpatt+1
    pat<-tabulate(intpatt,nbins=2^ncol(as.matrix(x)))
    pat
}
