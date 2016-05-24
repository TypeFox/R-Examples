`pair.prod` <-
function(X)
    {
     d<-dim(X)
     matrix(.C("pairprod", as.double(X),as.integer(d), res=double(choose(d[1],2)*d[2]),PACKAGE="ICSNP")$res,ncol=d[2],byrow=T)
    }
