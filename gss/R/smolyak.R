smolyak.quad <- ## Generate delayed Smolyak cubature
function(d, k) {
    size <- .C("size_smolyak",
               as.integer(d),
               as.integer(d+k),
               size=integer(1),
               PACKAGE="gss")$size
    z <- .C("quad_smolyak",
            as.integer(d),
            as.integer(d+k),
            pt=double(d*size),
            wt=as.double(1:size),
            PACKAGE="gss")
    list(pt=t(matrix(z$pt,d,size)),wt=z$wt)
}

smolyak.size <- ## Get the size of delayed Smolyak cubature
function(d, k) {
    .C("size_smolyak",
       as.integer(d),
       as.integer(d+k),
       size=integer(1),
       PACKAGE="gss")$size
}
