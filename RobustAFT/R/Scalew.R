"Scalew" <-
function(vo,rr,den,rhs,tol,maxit){
z  <- Nwtv(vo/4,rr,den,rhs,tol,maxit)$v
if (!is.na(z)) v <- z else v <- Fxdv(vo,rr,den,rhs,tol,300)$v; v}

