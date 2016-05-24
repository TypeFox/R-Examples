"F0w" <-
function(u,tol=0.0001,maxit=150) {# Fortran version of F0w.s
z  <- .Fortran("srf0w",u=as.double(u),tol=as.double(tol),maxit=as.integer(maxit),p=double(1))
z$p}

