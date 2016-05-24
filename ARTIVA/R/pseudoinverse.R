pseudoinverse <-
function (m, tol)
{
    msvd = fast.svd(m, tol)
    
    if (length(msvd$d) == 0)
    {
       return(
            array(0, dim(m)[2:1])
            )
    }
    else
    {
       return( 
            msvd$v %*% (1/msvd$d * t(msvd$u))
            )
     }    
}
