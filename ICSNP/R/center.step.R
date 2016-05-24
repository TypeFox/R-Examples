`center.step` <-
function(datas,A,maxiter,eps.center)                                            
    {
     spatial.median(datas%*%t(A),init=NULL,maxiter=maxiter,eps=eps.center,print.it=FALSE) %*% solve(A)
    }
