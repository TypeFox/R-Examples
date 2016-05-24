vdecomp <- function(x,v1,v2) {
    w1 <- v1 - project(v1,v2) ; w2 <- v2 - project(v2,v1) 
    p1 <- project(x,w1,type='v') ; p2 <- project(x,w2,type='v') 
    a  <- sign( dot(w1, p1) ) * vlength(p1) / vlength(w1)
    b  <- sign( dot(w2, p2) ) * vlength(p2) / vlength(w2)
    structure( list( coefficients=c(a,b), 
        x.orig = x, x.proj = a * v1 + b * v2 
        ))
}
