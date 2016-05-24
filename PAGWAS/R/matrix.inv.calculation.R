matrix.inv.calculation <-
function(V){
    cV=chol(V)
    V.inv=chol2inv(cV)
    rm(cV)
    V.inv
}
