
mubz<-function(category="copula",...) {
    do.call(paste("mubz.",category,sep=""),list(...))
}

v.star<-function(m,z=rep(1,m),...) {
    sum(z*Vectorize(function(i) mubz(b=(1:m)==i,
                                     z=z,
                                     ...)
                    )(1:m))
}

