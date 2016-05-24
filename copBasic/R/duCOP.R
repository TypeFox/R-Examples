"duCOP" <-
function(u, v, cop=NULL, para=NULL, ...) {
    #str(cop)
    return(u + v - cop(u,v, para=para, ...))
}
