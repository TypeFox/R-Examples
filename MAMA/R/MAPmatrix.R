MAPmatrix<-function(value.dis)
{
res<-ratio(value.dis)
unique.pat <- unique(res$X.string)
n.soft <- patternMatch(value.dis,unique.pat)
n.strong <- patternMatch.strong(value.dis,unique.pat)
unique.X <- patternmatrix(unique.pat,ncol(value.dis))
n.sig <- apply(unique.X,1,sum)
mat<-data.frame(unique.pat, n.soft, n.strong,n.sig)
return(mat)
}
