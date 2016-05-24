Ainv <-
function(GAB, x, tol=1e-12 )
  {
    ###  need something to make up for the lame-o
    ##   matlab code that does this h = G\x to get the inverse
    if(missing(tol)) tol = .Machine$double.eps
#    gsvd=svd(GAB,LINPACK=FALSE)
#    desolve=(gsvd$v%*%(diag(1/gsvd$d))%*%t(gsvd$u))%*%x
#    desolve = solve(t(GAB)%*%GAB ,  t(GAB) %*% x , tol=defaulttol)  
 ###    qg=base::qr(GAB,tol=defaulttol,LAPACK=TRUE)

    desolve= qr.solve(GAB, x, tol = tol)

    
    ### desolve=solve(qg,x)
    
    return(desolve)   
  }
