"surfuncCOP" <-
function(u, v, cop=NULL, para=NULL, ...) {
    return(1 - u - v + cop(u,v, para=para, ...))
}
