statdistr <-
function(tmat)
{
one<-round(matrix(1,nrow(tmat),ncol(tmat)))

I<-diag(nrow(tmat))
d<-rep(1,nrow(tmat))

t(d)%*%solve(I-tmat+one)
}
