# aus e1071
# adapted for ternary patterns
countpattern3<-function (x)
{
    intpatt<-apply(as.matrix(x),1,tern2int)
    intpatt<-intpatt+1
    pat<-tabulate(intpatt,nbins=3^ncol(as.matrix(x)))
    pat
}
