matrix.invroot.calculation <-
function(V){
    E=eigen(V)
    E.vector=E$vectors
    V.inv.root=E.vector%*%diag(1/sqrt(E$values))%*%t(E.vector)
    rm(E.vector,E)
    V.inv.root
}
